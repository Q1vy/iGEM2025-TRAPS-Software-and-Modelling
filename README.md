# iGEM2025-TRAPS-Software
Repository for Software and Modelling done for the TU-Dresden team for iGEM2025 competition.
General topic: RNA-binding Proteins 
- predicting binding affinities
- sequence comparison for target selection


Content description:
- Main target search: given a target sequence (of an RNA to bind), can find subsequences of a given length (binding site) that are as unique compared to reference sequences (transcriptom). Next goals: integrate gene annotations (name, description) and transcript frequency data (to weight the importance of off-target interactions)
- sequence data: annotated transcriptome, no concentration/frequency data yet. Wrong mCherry sequence, will be updated soon
- modeling: simple Boltzman statistics, no condensates or polymer physics, to predict the binding interactions of RNA binding proteins with their specific target, and off-target interactions with the transcriptome.
