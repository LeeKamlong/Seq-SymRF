File content description:
  disease name.txt--339 diseases name including in our dataset
  miRNA name.txt--917 miRNAs name including in our dataset
  MDA.mat-- a binary matrix that represents the known miRNA-disease association from HMDD3.0
  DiseaseSymScore.mat-- a matrix contain 322 symptoms and their score for 399 disease
  miRfeature.m-- a function that was calculated the feature for miRNA base on the secondary structure
  miRNASecondStructure.mat--917 miRNAs with their sequence and secondary structure
Software dependencies:
  Operating system: Windows
  Programming language: Matlab
  Other requirements: Matlab2018b or higher
  Memory: 8 GB RAM
How to use:
  (1)Download miRfeature.m to the root directory of Matlab
  (2)Run the Seq-SymRF.m to construct the RF model and test the performance
  (3)Save the RF model
How to predict the miRNA or disease which I want:
  For miRNAs:
    (1)You should obtain the pre-miRNA sequence for miRbase (eg: miR-21-5p and miR-21-3p correspond to the same pre-miRNA sequence hsa-mir-21)
    (2)Use the wed sever RNAfold (http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi) to predict the secondary sequence
    (3)Run the miRfeature.m to calculate the feature for miRNA
    (4)Combine the miRNA feature with disease feature
    (5)Prediction
  For diseases:
    (1）Obtain the symptom score from the disease-symptom network
    (2) Combine the disease feature with miRNA feature
    (3)Prediction
 %% If you have any questions, please contact me zhanchao8052@gdpu.edu.cn (Prof.Li) and ceszxy@mail.sysu.edu.cn (Prof.Zhou).
