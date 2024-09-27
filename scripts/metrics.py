# pip install pandas

import pandas as pd


df = pd.read_table(snakemake.input[0])      # important to use the snakemake structure. Check double brackets 
outFile = open(snakemake.output[0], 'w')    # snakemake structure and double brackets. 


#print(df.head(5))
#print(df.info()) ### infos about dataframe

#print(len(df))


# Variable definitions
samplesOccurrence = df.groupby('Sample').size()
classOccurrence = df.groupby('Class').size()
timepointsOccurrence = df.groupby('Timepoint').size()
animalsOccurrence = df.groupby('Animal').size()

# Writing in a new file
outFile.write("Número de amostras: " + "\n" + str(samplesOccurrence))
outFile.write("\n\n" + "Ocorrência por classe de amostra: " + "\n" + str(classOccurrence))
outFile.write("\n\n" + "Ocorrência por time points: " + "\n" + str(timepointsOccurrence))
outFile.write("\n\n" + "Número de animais: " + "\n" + str(animalsOccurrence))


outFile.close()
