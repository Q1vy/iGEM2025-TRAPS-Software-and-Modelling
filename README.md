# iGEM2025-TRAPS-Software
Repository for Software and Modelling done for the TU-Dresden team for iGEM2025 competition.
General topic: RNA-binding Proteins 
- predicting binding affinities
- sequence comparison for target selection


Content description:
- Main target search: given a target sequence (of an RNA to bind), can find subsequences of a given length (binding site) that are as unique compared to reference sequences (transcriptom).
  - Currently done with a preliminary mCherry sequence as target (will be updates soon) tested against the transcriptome of Sacch. Cerevisae. Total rununtime < 10min.
  - Next goals: integrate gene annotations (name, description) and transcript frequency data (to weight the importance of off-target interactions)
- sequence data: annotated transcriptome, no concentration/frequency data yet. Wrong mCherry sequence, will be updated soon
- modeling: simple Boltzman statistics, to predict the binding interactions of RNA binding proteins with their specific target, and off-target interactions with the transcriptome. Possible next steop: include aspects of condensates or polymer physics
