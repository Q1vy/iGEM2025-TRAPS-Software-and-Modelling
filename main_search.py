
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import pandas as pd
import warnings
import pathlib




@njit  # just in time compilation speeds up sequence comparison, also see seq_compare_runtime file
def sc_jit(seq:bytes|str, ref:bytes|str):
    '''A fast Sequence Comparison of a short query sequence (seq) against a larger reference (ref).
    Returns the number of matches for every position the sort sequence could bind to the larger one. Consideres only full overlap binding, i.e. no overhangs '''

    n = len(ref) - len(seq) + 1
    if n<1: print("sc_jit: wrong argument order")
    cv = np.zeros(n, dtype=np.uint16) 
    for i in range(n):
        matches = 0
        for j in range(len(seq)):
            if ref[i + j] == seq[j]:
                matches += 1
        cv[i] = matches
    return cv

def random_seq(size:int):
    '''returns a random ACGT sequence of length = size'''
    return ''.join(np.array(["A", "C", "G", "T"])[np.random.randint(low = 0, high = 4, size = size)])

sc_jit(random_seq(10), random_seq(10));    # call it once so that numba compiles it




def query_sc(query:pd.Series, reference_dataset:pd.DataFrame)->pd.DataFrame:
    '''searches for (partial) matches between query sequence and reference dataset
    query: sequence to be found (pd.Series object)
    reference_dataset: pd.DataFrame of sequences to test against
    Both need at least a SEQ_NAME and a SEQUENCE column
    Returns copy of reference_dataset with new columns containing the results
    '''
    #-------------------------------------------------#
    #              input arguments check              #  -> checks if the parameters have expected datatype to prevent errors later on
    #-------------------------------------------------#

    # check data types
    if not type(reference_dataset) == pd.DataFrame: raise TypeError("unsupported type for 'reference_dataset', use pandas.DataFrame")
    if not ("SEQ_NAME" in reference_dataset.columns and "SEQUENCE" in reference_dataset.columns):
        raise TypeError('reference_dataset does not contain required columns: SEQ_NAME, SEQUENCE')
    if type(query) == pd.Series:  pass # everything is fine
    elif type(query) == str:           # convert to pandas series
        warnings.warn("Warning, deprecated input type for query: string. Please use pandas.Series!")
        query = pd.Series({"SEQ_NAME": "unnamed", "SEQUENCE": query})
    elif type(query) == pd.DataFrame and len(query) == 1:
        warnings.warn("Warning, query had to be converted from pandas.DataFrame to pandas.Series!")
        query = query.iloc[0]
    else: raise TypeError("unsupported type for 'query', use pandas.Series")

    # check how many different nucleotides occure in the provided sequences
    nts = set(reference_dataset.SEQUENCE.sample(100).sum()) | set(reference_dataset.SEQUENCE.iloc[:10].sum()) | set(reference_dataset.SEQUENCE.iloc[-10:].sum()) | set(query.SEQUENCE)
    if len(nts) != 4: 
        raise ValueError(f"Warning, expected 4 types of nucleoties in the input, found (at least) {len(nts)}")

    #-------------------------------------------------#
    #  Search for binding sites in reference dataset  #
    #-------------------------------------------------# 
    sc_results = []                                                          # initialize a dataset to store the results of comparing query sequence with every reference i.e. every transcript in the transcriptome
    for i, ref in reference_dataset.iterrows():                              # iterate over all reference sequences; and the target itself as control 

        # check how many different nucleotides occure in the provided sequences
        nts = len(set(ref.SEQUENCE) | set(query.SEQUENCE) )
        if nts != 4: raise ValueError(f"Warning, expected 4 types of nucleoties in the input, found (at least) {len(nts)}")


        # number of mismatches per position  = length - number of matches. I dunno why but its faster in bytes
        if type(query.SEQUENCE) == str:     mm = len(query.SEQUENCE)-sc_jit(query.SEQUENCE.encode(), ref.SEQUENCE.encode())    
        elif type(query.SEQUENCE) == bytes: mm = len(query.SEQUENCE)-sc_jit(query.SEQUENCE, ref.SEQUENCE)      

        sc_row = {}                                                          # create dictionary to summarize the possible binding sites of the query sequence against one reference sequence (here a transcript). sc = sequence compare
        sc_row["identical match"]  = (mm==0).sum()                           # number of idential matches
        sc_row["1nt mismatch"]     = (mm==1).sum()                           # number of sites with 1 mismatch
        sc_row["2nt mismatch"]     = (mm==2).sum()                           # number of sites with 2 mismatch
        sc_row["boltzmann factor"] = np.sum(np.exp(-2*mm.astype(float)))     # exp(energy) with energy ~ 2kt * mm is a very basic prediction of how likely the query sequence would bind to the reference. The wildtype pumilios have ~ 18 hydrogen bonds for 8nt i.e. 2.25 hydgrogen bonds = 2.25kT per bp, however the pumby appears to only have ~2 H-bonds per bp, and the binding statistics suggest closer to 1kT -> lets go with 1.5 for now, in the future a more accurate model might be implemented. Currently the boltzmann_factor for a perfect match is 1, if therefore this predicted boltzmann_factor is smaller than 1 than it is more likly to bind to the specific RNA, while if it is larger than one it is more likely to bind to the background transcriptom. The float conversion is there because pandas does not like float16 only float32 and larger, for some reason without this explicit astype it would be float16 
        sc_results.append(sc_row)                                            # add the sequence comparison summary for this reference (here transcript) to the dataset
    sc_results = pd.concat([reference_dataset, pd.DataFrame(sc_results, reference_dataset.index)], axis = 1)   # convert to pandas dataframe and add original columns 

    return sc_results




def round_format(x:np.number, digits: int = 3):
    '''round numbers to "digits" significant digits. i.e. 1234 with digits = 2 -> 1200, and 0.1234 -> 0.12'''
    if x == 0: return 0
    mag = -int(np.floor(np.log10(x))-digits+1)
    return np.round(x, mag).astype([int, float][int(mag > 0)])


def sc_summarize(query:pd.Series, references_scr:pd.DataFrame, target_scr:pd.Series, use_weights = True, VIR_max = 10, r2t = True, plot = False) -> pd.Series:
    ''' function summarizing the result table from query_sc'''
    #-------------------------------------------------#
    #              input arguments check              #  -> checks if the parameters have expected datatype to prevent errors later on
    #-------------------------------------------------#

    if type(references_scr) != pd.DataFrame: raise TypeError("unsupported type for 'references_scr', use pandas.DataFrame")

    if use_weights:  # check if the input data is okay to use weighted evalutation
        if "WEIGHT" in references_scr.columns:
            invalid = references_scr.WEIGHT.isin([np.inf, -np.inf, np.nan])
            if invalid.any():  warnings.warn(f"{invalid.sum()} entries with nan or inf WEIGHT in references_scr. The script might not behave as intended.")
        else:   raise TypeError("references_scr: No WEIGHT column found, to switch to unweighted mode use use_weights = False")

    if r2t:
        if type(target_scr) == pd.Series:  pass
        elif type(target_scr) == pd.DataFrame:  
            warnings.warn("Warning, query had to be converted from pandas.DataFrame to pandas.Series!")
            target_scr = target_scr.iloc[0]
        else: raise TypeError("unsupported type for 'target_scr', use pandas.Series")

        if use_weights:  # check if the input data is okay to use weighted evalutation
            if "WEIGHT" in target_scr:
                if target_scr.WEIGHT in [np.inf, -np.inf, np.nan]: warnings.warn(f"target_scr has nan or inf WEIGHT. The script might not behave as intended.")
            else:   raise TypeError("target_scr: No WEIGHT property found, to switch to unweighted mode use use_weights = False")



    #-------------------------------------------------------#
    # summarize interaction of query with reference dataset #   ->  now that we have a huge and detailed list of the interaction with every of the thousands of reference sequences (i.e. transcripts), we need to summarize again to only a few metrics telling how much this query sequence interacts with the references overall
    #-------------------------------------------------------#

    sc_summary = query.to_dict()
    for i, metric in enumerate(["identical match", "1nt mismatch", "2nt mismatch", "boltzmann factor"]): # iterate over all the metrics that should be summarized

        ref_metric = references_scr[metric]         # results of comparing the query i.e. candidate sequence against the Reference dataset. metrics from references -> ref_metric
        if r2t: trg_metric = target_scr[metric]             # results of comparing the query i.e. candidate sequence against the Target sequence.   metrics from target   -> trg_metric
        if use_weights:   # weighted metric -> wetric
            ref_wetric = references_scr[metric] * references_scr["WEIGHT"]
            if r2t: trg_wetric = target_scr[metric] * target_scr["WEIGHT"]

        if "oltzman" in metric and r2t:               # for the Boltzmann factor it is compared to the Boltzmann factor of the query sequence --> relative comparison of query binding to the target vs transcripts
            insert = " r2t"; norm = trg_metric        # = "relative to target". used for column name
            if use_weights: norw = trg_wetric         # norw = normalization for weighted data
        else:  insert = "";   norm = 1;   norw = 1    # the mismatch count data is not re-normalized

        # summarize the sequence comparison results: sum = total number of binding sites, while max = maximum number of binding sites per reference sequence i.e. transcript
        sc_summary[ metric + f" max" + insert] = ref_metric.max() / norm        # What was the highest number of off-target binding sites on a single reference i.e. transcript
        sc_summary[ metric + f" sum" + insert] = ref_metric.sum() / norm        # What was the  total  number of off-target binding sites on the whole reference dataset i.e. transcriptome
        VIR_selec = ((ref_metric > sorted(ref_metric)[-min(VIR_max, len(ref_metric))]) | (ref_metric == ref_metric.max()) * (ref_metric != ref_metric.min())) ## (ref_metric > sorted(ref_metric)[-10]) returns the 10th largest value, also ever value equal to the largest one is included. Overrulling everything else is the exclusion of the min value
        VIR_names = references_scr.iloc[ref_metric[VIR_selec].sort_values(ascending = False).index].SEQ_NAME  # Store Very Important References i.e. those with strong offtarget binding sites
        sc_summary[ metric + " VIR"] = VIR_names.tolist()
        if use_weights:
            sc_summary[ metric + f" W.max" + insert] = ref_wetric.max() / norw        # What was the highest number of off-target binding sites on a single reference i.e. transcript
            sc_summary[ metric + f" W.sum" + insert] = ref_wetric.sum() / norw        # What was the  total  number of off-target binding sites on the whole reference dataset i.e. transcriptome
            VIR_selec = ((ref_wetric > sorted(ref_wetric)[-min(VIR_max, len(ref_wetric))]) | (ref_wetric == ref_wetric.max()) * (ref_wetric != ref_wetric.min())) ## (ref_metric > sorted(ref_metric)[-10]) returns the 10th largest value, also ever value equal to the largest one is included. Overrulling everything else is the exclusion of the min value
            VIR_names = references_scr.iloc[ref_wetric[VIR_selec].sort_values(ascending = False).index].SEQ_NAME  # Store Very Important References i.e. those with strong offtarget binding sites
            sc_summary[ metric + " W.VIR"] = VIR_names.tolist()


        #------------------------------------------------------#
        #                 Plotting the results                 #   
        #------------------------------------------------------#


        # still within the same for loop iterating over every metric 
        if plot:


            def round_format(x, digits = 3):
                '''round numbers to "digits" significant digits. i.e. 1234 with digits = 2 -> 1200, and 0.1234 -> 0.12'''
                if x == 0: return 0
                mag = -int(np.floor(np.log10(x))-digits+1)
                return np.round(x, mag).astype([int, float][int(mag > 0)])

            xmax = len(ref_metric) # number of references
            sy = 1  # weighted data is shown in the same plot, but needs a seperate, rescaled y-axis. sy is this rescaling factor
            if use_weights: 
                if r2t: _data = pd.concat([references_scr, pd.DataFrame([target_scr ])], axis = 0)
                else:   _data = references_scr
                if _data[metric].max() != 0 and (_data[metric]*_data.WEIGHT).max() != 0 and (_data[metric]*_data.WEIGHT).min() >= 0: # avoid divide by zero in case of no data 
                    sy = _data[metric].max() / (_data[metric]*_data.WEIGHT).max() # internal rescaling factor for the split y axis
                if r2t:
                    sy /= norm/norw

            fig, ax = plt.subplots(1, 1, figsize = (12, 3 + .5*int(use_weights)))
            ax.plot(        [0, xmax],                [0, 0], zorder = 10, color = "black", lw = 1) # plot black zero line for references
            if r2t: ax.plot([-.032*xmax, -.028*xmax], [0, 0], zorder = 10, color = "black", lw = 1) # plot black zero line for target
            ax.fill_between(        np.arange(xmax),           ref_metric/norm,              ref_metric*0, color = "tab:blue", ec = "tab:blue", lw = .5, label = "unweighted")    # plot unweighted results of references
            if r2t: ax.fill_between([-.032*xmax, -.028*xmax], [trg_metric/norm, trg_metric/norm], [0,0],   color = "tab:blue", ec = "tab:blue", lw = .5)                          # plot unweighted results of target

            if use_weights:
                ax.fill_between(        np.arange(xmax),           -ref_wetric*sy/norw,                   ref_wetric*0,  color = "tab:orange", ec = "tab:orange", lw = .5, label = "weighted")           # plot weighted results of references
                if r2t: ax.fill_between([-.032*xmax, -.028*xmax], [-trg_wetric*sy/norw, -trg_wetric*sy/norw], [0,0],     color = "tab:orange", ec = "tab:orange", lw = .5)                               # plot weighted results of target
                ax.legend(loc = "upper right")

            message = "" # for figure titles
            if "match" in metric: # for the (mis) match count data
                ax.set_ylabel("number of binding sites")
                if ref_metric.max() > 0:  # i.e. if there are mismatches to plot.
                    message = f"  »  max {ref_metric.max()}x on {(ref_metric == ref_metric.max()).sum()} reference{'s'*int((ref_metric == ref_metric.max()).sum() != 1)}"
            elif "oltzma" in metric: # for the boltzmann factors
                if use_weights:  message =  f" »  w.sum≈{round_format(ref_wetric.sum()/norw)} " + insert
                else:            message =  f" »  sum≈{  round_format(ref_metric.sum()/norm)} " + insert
                if r2t: ax.set_ylabel("relative binding frequency")
                else:   ax.set_ylabel("binding frequency [au]")

            ax.set_title(metric + message)

            # x-axis formatting
            ax.set_xlabel("reference sequence name")
            if r2t: 
                ax.set_xticks(np.concatenate([[-0.03*xmax], VIR_names.index]), labels = np.concatenate([["target"], VIR_names.tolist()]), fontsize = 8, rotation = 45,  ha='right')   
                ax.set_xlim(-0.05*xmax, 1.01*xmax)

            else:  
                ax.set_xticks(VIR_names.index, labels =  VIR_names.tolist(), fontsize = 8, rotation = 45,  ha='right')
                ax.set_xlim(-0.01*xmax, 1.01*xmax)

            # split y axis formatting
            yt = ax.get_yticks().copy()
            ax.set_yticks(yt[yt>=0])
            sec_y = ax.secondary_yaxis("left", functions = (lambda y: -y/sy, lambda y: -y*sy))
            plt.tight_layout() # need to be called before getting ytick values as this triggers matplotlib to calculate the layout and the ticks correctly.
            yt2 = sec_y.get_yticks()
            sec_y.set_yticks(yt2[yt2>0])

            plt.show()

    return pd.Series(sc_summary)




def query_sc_summary(query:pd.Series, reference_dataset:pd.DataFrame, target:pd.Series = None, use_weights: bool = True, VIR_max: int = 10, r2t: bool = True, plot: bool = False):
    '''searches for (partial) matches between query sequence and reference dataset to estimate binding affinities. If query is supposed to bind to a target, 
    > query: sequence to be found (pd.Series object)
    > reference_dataset: pd.DataFrame of sequences to test against i.e. transcriptome
    > target (optional): a single reference that the query is supposed to bind to. formatted like query.
    Both need at least a SEQ_NAME and a SEQUENCE column
    '''

    ### Input handling ### -> most of it is done by "query_sc"
    if r2t:
        if type(target) == pd.Series: sc_refs = pd.concat([pd.DataFrame([target]), reference_dataset], axis = 0)  # if we use the target, then add it to the references used in the sequence comparison
        else: raise TypeError("target: expected pd.Series")
    else: sc_refs = reference_dataset

    # run outsourced sequence comparison
    sc_results = query_sc(query=query, reference_dataset=sc_refs) # reference_dataset in the "query_sc" is the dataset of references to test against, which can include the target, in the big picture reference_dataset is only references like transcripts

    # retrieve target and references results
    if r2t:
        references_scr = sc_results.loc[~(sc_results.SEQ_NAME == target.SEQ_NAME)].reset_index(drop = True)  # results from the comparison with the reference dataset (sequences to avoid)
        _target_scr =    sc_results.loc[  sc_results.SEQ_NAME == target.SEQ_NAME ].reset_index(drop = True)  # results from the comparison with the target sequence
        if len(_target_scr) != 1: raise IndexError(f"unexpected result, there should be one entry, found {len(_target_scr)}")      # there should only be one entry for the target
        else: target_scr = _target_scr.iloc[0]
    else: 
        references_scr = sc_results
        target_scr = None

    # run outsourced summarizer
    sc_summary = sc_summarize(query, references_scr=references_scr, target_scr=target_scr, use_weights=use_weights, VIR_max=VIR_max, r2t=r2t, plot=plot)

    return sc_summary



import scipy.stats as stats

def fig_dim(n):
    r = n/(1.25*n**.5)
    y = np.ceil(r * (r**2/(r**2+5)))
    x = np.ceil(n/y)
    ncol = x
    nrow = y
    for o in np.arange(-0.33*x, 0.5*x).round():
        _x = x+o
        _y = math.ceil(n/_x)
        if (_x*_y-n < nrow*ncol-n):
            ncol = x
            nrow = y
    return int(min(nrow, ncol)), int(max(ncol, nrow))

def df_PCA(df, show_cols = None, exclude_cols = ["POSITION"]):
    '''plots results of all queries in a PCA scatterplot. cols = DAtaFrame columns to plot. All numeric columns are used for PCA.
    show_cols = columns to plot. Default = None = all columns
    exclude_cols = columns to exclude from PCA. All non-numeric columns are ignored either way.'''
    df = df.copy()
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    exclude_cols = [col for col in df.columns if np.any([exc in col for exc in exclude_cols])]
    data = StandardScaler().fit_transform(df.select_dtypes(np.number).drop(exclude_cols, axis = 1))
    PCA1, PCA2 = PCA(n_components=2).fit_transform(data).T
    dfn = df.select_dtypes(np.number)

    if type(show_cols) == type(None):
        columns = dfn.columns
    else:
        columns = [col for col in show_cols if col in dfn.columns]

    nrow, ncol = fig_dim(len(columns))

    fig, axes = plt.subplots(nrow, ncol, figsize = (2.5*ncol, 2*nrow))
    for col, ax in zip(columns, np.array([axes]).flatten()):
        points = ax.scatter(PCA1, PCA2, c = dfn[col], s = 5, alpha = 1, cmap = ["Spectral_r", "Spectral"][int(col == "POSITION")]) # Position has the cmap the other way around to be consistent with the table, while all other numbers use _r so that high = bad numbers = red similar to table again
        ax.set_title(col)
        ax.set_yticks([])
        ax.set_xticks([])
        plt.colorbar(points, ax = ax)
    plt.tight_layout()
    plt.show()

def df_map(df, columns = ["identical match", "1nt mismatch", "2nt mismatch", "boltzmann factor"]):
    '''Visualize the number of offtarget binding sites and boltzmann factors of all queries from the target sequence'''
    df = df.copy()
    for col in columns:
        _cols = df.columns[np.array([col in _c for _c in df.columns])]
        col1 = _cols[["sum" in _c for _c in _cols]][0]
        col2 = _cols[["max" in _c for _c in _cols]][0]
        print(col1)
        _min = 0
        if "W." in col1: _min = 0
        sel = df[df[col2]==_min]

        plt.figure(figsize = (20, 3))
        sy = np.round(0.66 * df[col1].mean() / df[col2].mean(), 0).astype(int) # rescaling for the lower half of the plot, only visual
        plt.bar(df.POSITION,  df[col1].values, label = "total sum in transcriptome", width = 1, color = "tab:cyan")
        plt.bar(df.POSITION, -df[col2].values*sy,  label = "max per transcript", width = 1, color = "tab:purple")
        if len(sel) > 0 and df[col2].min() == 0:
            plt.bar(sel.POSITION, .04*df[col1].max(), bottom = -.02*df[col2].max()*sy, label = f"no {col} off-target binding sites in any transcript", color = "black", width = 1)
        plt.xlabel("target candidates")
        yl = plt.yticks()[0]
        yl[yl<0] /= sy
        yl = np.unique(yl.round().astype(int))
        yt = yl.copy()
        yt[yt<0] *= sy
        plt.yticks(yt, np.abs(yl))
        plt.legend();

        if not  "oltzman" in col:
            plt.ylabel(f"number of {col} \n off-target binding sites");

            if len(sel) == 0:
                plt.title(f"no candidate without {col} off-target binding sites found")
            elif df[col2].min() == 0:
                plt.title(f"{len(sel)} candidates without {col} off-target binding sites found")
            else:
                plt.title(f"{len(sel)} candidates with max {round(df[col2].min())}x {col} off-target binding sites found")
        else: 
            plt.title("relative boltzmann factors normalised to BF of target sequence")
            plt.ylabel(f"{col} of \n off-target binding sites");


        plt.show()

def _col_weights(col):
    col_name = col.name
    if "W." in col_name: m = 2
    else: m = 1
    if "identical" in col_name or "boltzman" in col_name:    return m*1
    elif "1nt" in col_name:                                  return m*1/2
    elif "2nt" in col_name:                                  return m*1/4
    else:                                                    return m*1/8

def nice_format(var, digits = 3):

        try:
            x = float(var)
            if var == 0: return int(0)
            if np.isnan(var): return "---"
            mag = -int(np.floor(np.log10(x))-digits+1)
            num = np.round(x, mag)
            if mag > 0:
                if int(var) == var: return int(var)
                return (f"{num:.{mag}f}")
            else: 
                return (f"{int(num)}")    
        except: return var

def df_table(DF, select_by_max = "selected", ignore_columns_like = ["VIR", "WEIGHT"]):

    col_names = list(DF.columns)
    levels = col_names.copy()
    for i, col in enumerate(col_names):
        candidate = col
        found = False
        while True:
            d =  max([candidate[:-1].rfind(d) for d in [" "]])
            if d == -1: break; # no more splitting possible
            candidate = col[:d+1]
            if sum([s.startswith(candidate) for s in col_names]) > 1: found = True; break; # found common column name substring
        if found:
            k = col.find(candidate) + len(candidate)
            levels[i] = (col[:k], col[k:])
        else: levels[i] = (col, "")
    levels  

    df = DF.copy().sort_index()
    if type(select_by_max) == str: 
        if select_by_max in DF.columns:
            df = DF[DF[select_by_max] == DF[select_by_max].max()].copy().sort_index()

    median = DF.select_dtypes("number").median()
    median[[col for col in df.columns if col.startswith("POSITION")]] = np.nan 
    _df = pd.concat([df, pd.DataFrame([median], index = [None])], axis = 0, join = "outer")
    _df.columns = pd.MultiIndex.from_tuples(levels) 

    position_cols = [col for col in _df.columns if col[0].startswith("POSITION")]

    header0  = {'selector': 'th.col_heading.level0', 'props': 'font-size: 1.5em;'}
    header1  = {'selector': 'th.col_heading.level1', 'props': 'border-bottom: 5px solid white'}
    footer  = {'selector': '.data', 'props': 'border-top: 2px solid white'}
    general = {'selector': 'th, td',                     'props': 'text-align: center;'}
    #_df = _df.reset_index(names = ["position"])

    _s = _df.style.background_gradient(axis=0, cmap = "Spectral", subset=position_cols, vmin = 0, vmax = DF.filter(like = "POSITION").max().max() ).hide(axis='index')

    for col in _df.drop(position_cols, axis = 1).select_dtypes("number").columns:
        _s = _s.background_gradient(axis=0, cmap = "coolwarm", subset = [col], vmax = (DF[col[0] + col[1]].median()*2 - DF[col[0] + col[1]].min()) + 1e-12, vmin = (1*DF[col[0] + col[1]].min()) - 1e-12 )
    _s = _s.hide([col for col in _df.columns if np.any([pattern in col[1] for pattern in ignore_columns_like]) or np.any([pattern in col[0] for pattern in [select_by_max]+ignore_columns_like])], axis = 1,)
    _s = _s.format(nice_format).format((lambda x: f"{x:.0f}" if not np.isnan(x)  else "DF.Median"), subset = position_cols).set_table_styles([header0, header1, general, footer])
    _s = _s.set_table_styles({head: [{'selector': 'th, td', 'props': 'border-left: 5px solid white'}] for head in _df.columns if (head[1].startswith("max")) or (head[1] == "")}, overwrite=False, axis=0)
    _s = _s.set_table_styles({pc:[{'selector': 'td', 'props': "font-weight: bold"}] for pc in position_cols},  overwrite=False, axis=0)
    _s = _s.set_table_styles({("SEQUENCE", ""):[{'selector': 'td', 'props': "font-family: monospace"}]},  overwrite=False, axis=0)
    display(_s)

def flatten(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))  # Recursively flatten
        else:
            flat_list.append(item)
    return flat_list

def str_multisplit(string, substring_list):
    for k, delimiter in enumerate(substring_list):
        if k == 0: split = string.split(delimiter)
        else: split = flatten(s.split(delimiter) for s in split)
    return split

def str_alt_split(string, substring):
    '''splits at the end of the substring without removing the substring, but trimming the split results'''
    i = string.find(substring) + len(substring)
    if i == -1:
        return string
    else: return [string[:i].strip(), string[i:].strip()]

import scipy.stats as stats
import math

def fig_dim(n):
    r = n/(1.25*n**.5)
    y = np.ceil(r * (r**2/(r**2+5)))
    x = np.ceil(n/y)
    ncol = x
    nrow = y
    for o in np.arange(-0.33*x, 0.5*x).round():
        _x = x+o
        _y = math.ceil(n/_x)
        if (_x*_y-n < nrow*ncol-n):
            ncol = x
            nrow = y
    return int(min(nrow, ncol)), int(max(ncol, nrow))

def df_PCA(df, show_cols = None, exclude_cols = ["POSITION"]):
    '''plots results of all queries in a PCA scatterplot. cols = DAtaFrame columns to plot. All numeric columns are used for PCA.
    show_cols = columns to plot. Default = None = all columns
    exclude_cols = columns to exclude from PCA. All non-numeric columns are ignored either way.'''
    df = df.copy()
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    exclude_cols = [col for col in df.columns if np.any([exc in col for exc in exclude_cols])]
    data = StandardScaler().fit_transform(df.select_dtypes(np.number).drop(exclude_cols, axis = 1))
    PCA1, PCA2 = PCA(n_components=2).fit_transform(data).T
    dfn = df.select_dtypes(np.number)

    if type(show_cols) == type(None):
        columns = dfn.columns
    else:
        columns = [col for col in show_cols if col in dfn.columns]

    nrow, ncol = fig_dim(len(columns))

    fig, axes = plt.subplots(nrow, ncol, figsize = (2.5*ncol, 2*nrow))
    for col, ax in zip(columns, np.array([axes]).flatten()):
        points = ax.scatter(PCA1, PCA2, c = dfn[col], s = 5, alpha = 1, cmap = ["Spectral_r", "Spectral"][int(col == "POSITION")]) # Position has the cmap the other way around to be consistent with the table, while all other numbers use _r so that high = bad numbers = red similar to table again
        ax.set_title(col)
        ax.set_yticks([])
        ax.set_xticks([])
        plt.colorbar(points, ax = ax)
    plt.tight_layout()
    plt.show()

def df_map(df, columns = ["identical match", "1nt mismatch", "2nt mismatch", "boltzmann factor"]):
    '''Visualize the number of offtarget binding sites and boltzmann factors of all queries from the target sequence'''
    df = df.copy()
    for col in columns:
        _cols = df.columns[np.array([col in _c for _c in df.columns])]
        col1 = _cols[["sum" in _c for _c in _cols]][0]
        col2 = _cols[["max" in _c for _c in _cols]][0]
        print(col1)
        _min = 0
        if "W." in col1: _min = 0
        sel = df[df[col2]==_min]

        plt.figure(figsize = (20, 3))
        sy = np.round(0.66 * df[col1].mean() / df[col2].mean(), 0).astype(int) # rescaling for the lower half of the plot, only visual
        plt.bar(df.POSITION,  df[col1].values, label = "total sum in transcriptome", width = 1, color = "tab:cyan")
        plt.bar(df.POSITION, -df[col2].values*sy,  label = "max per transcript", width = 1, color = "tab:purple")
        if len(sel) > 0 and df[col2].min() == 0:
            plt.bar(sel.POSITION, .04*df[col1].max(), bottom = -.02*df[col2].max()*sy, label = f"no {col} off-target binding sites in any transcript", color = "black", width = 1)
        plt.xlabel("target candidates")
        yl = plt.yticks()[0]
        yl[yl<0] /= sy
        yl = np.unique(yl.round().astype(int))
        yt = yl.copy()
        yt[yt<0] *= sy
        plt.yticks(yt, np.abs(yl))
        plt.legend();

        if not  "oltzman" in col:
            plt.ylabel(f"number of {col} \n off-target binding sites");

            if len(sel) == 0:
                plt.title(f"no candidate without {col} off-target binding sites found")
            elif df[col2].min() == 0:
                plt.title(f"{len(sel)} candidates without {col} off-target binding sites found")
            else:
                plt.title(f"{len(sel)} candidates with max {round(df[col2].min())}x {col} off-target binding sites found")
        else: 
            plt.title("relative boltzmann factors normalised to BF of target sequence")
            plt.ylabel(f"{col} of \n off-target binding sites");


        plt.show()

def _col_weights(col):
    col_name = col.name
    if "W." in col_name: m = 2
    else: m = 1
    if "identical" in col_name or "boltzman" in col_name:    return m*1
    elif "1nt" in col_name:                                  return m*1/2
    elif "2nt" in col_name:                                  return m*1/4
    else:                                                    return m*1/8

def nice_format(var, digits = 3):

        try:
            x = float(var)
            if var == 0: return int(0)
            if np.isnan(var): return "---"
            mag = -int(np.floor(np.log10(x))-digits+1)
            num = np.round(x, mag)
            if mag > 0:
                if int(var) == var: return int(var)
                return (f"{num:.{mag}f}")
            else: 
                return (f"{int(num)}")    
        except: return var

def df_table(DF, select_by_max = "selected", ignore_columns_like = ["VIR", "WEIGHT"]):

    col_names = list(DF.columns)
    levels = col_names.copy()
    for i, col in enumerate(col_names):
        candidate = col
        found = False
        while True:
            d =  max([candidate[:-1].rfind(d) for d in [" "]])
            if d == -1: break; # no more splitting possible
            candidate = col[:d+1]
            if sum([s.startswith(candidate) for s in col_names]) > 1: found = True; break; # found common column name substring
        if found:
            k = col.find(candidate) + len(candidate)
            levels[i] = (col[:k], col[k:])
        else: levels[i] = (col, "")
    levels  

    df = DF.copy().sort_index()
    if type(select_by_max) == str: 
        if select_by_max in DF.columns:
            df = DF[DF[select_by_max] == DF[select_by_max].max()].copy().sort_index()

    median = DF.select_dtypes("number").median()
    median[[col for col in df.columns if col.startswith("POSITION")]] = np.nan 
    _df = pd.concat([df, pd.DataFrame([median], index = [None])], axis = 0, join = "outer")
    _df.columns = pd.MultiIndex.from_tuples(levels) 

    position_cols = [col for col in _df.columns if col[0].startswith("POSITION")]

    header0  = {'selector': 'th.col_heading.level0', 'props': 'font-size: 1.5em;'}
    header1  = {'selector': 'th.col_heading.level1', 'props': 'border-bottom: 5px solid white'}
    footer  = {'selector': '.data', 'props': 'border-top: 2px solid white'}
    general = {'selector': 'th, td',                     'props': 'text-align: center;'}
    #_df = _df.reset_index(names = ["position"])

    _s = _df.style.background_gradient(axis=0, cmap = "Spectral", subset=position_cols, vmin = 0, vmax = DF.filter(like = "POSITION").max().max() ).hide(axis='index')

    for col in _df.drop(position_cols, axis = 1).select_dtypes("number").columns:
        _s = _s.background_gradient(axis=0, cmap = "coolwarm", subset = [col], vmax = (DF[col[0] + col[1]].median()*2 - DF[col[0] + col[1]].min()) + 1e-12, vmin = (1*DF[col[0] + col[1]].min()) - 1e-12 )
    _s = _s.hide([col for col in _df.columns if np.any([pattern in col[1] for pattern in ignore_columns_like]) or np.any([pattern in col[0] for pattern in [select_by_max]+ignore_columns_like])], axis = 1,)
    _s = _s.format(nice_format).format((lambda x: f"{x:.0f}" if not np.isnan(x)  else "DF.Median"), subset = position_cols).set_table_styles([header0, header1, general, footer])
    _s = _s.set_table_styles({head: [{'selector': 'th, td', 'props': 'border-left: 5px solid white'}] for head in _df.columns if (head[1].startswith("max")) or (head[1] == "")}, overwrite=False, axis=0)
    _s = _s.set_table_styles({pc:[{'selector': 'td', 'props': "font-weight: bold"}] for pc in position_cols},  overwrite=False, axis=0)
    _s = _s.set_table_styles({("SEQUENCE", ""):[{'selector': 'td', 'props': "font-family: monospace"}]},  overwrite=False, axis=0)
    display(_s)

def flatten(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))  # Recursively flatten
        else:
            flat_list.append(item)
    return flat_list

def str_multisplit(string, substring_list):
    for k, delimiter in enumerate(substring_list):
        if k == 0: split = string.split(delimiter)
        else: split = flatten(s.split(delimiter) for s in split)
    return split

def str_alt_split(string, substring):
    '''splits at the end of the substring without removing the substring, but trimming the split results'''
    i = string.find(substring) + len(substring)
    if i == -1:
        return string
    else: return [string[:i].strip(), string[i:].strip()]

