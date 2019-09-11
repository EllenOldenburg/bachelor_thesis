import csv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans
from sklearn import preprocessing



dataset_txt = "colombos_ecoli_exprdata_20151029.txt"
dataset_csv = "colombos_ecoli_exprdata_20151029.csv"
exp_dataset_txt = "colombos_ecoli_refannot_20151029.txt"
expid_dataset_txt = "colombos_ecoli_exprdata_20151029_2.txt"
expid_dataset_testannot_txt = "colombos_ecoli_testannot_20151029.txt"
expid_dataset_testannot_growth_txt = "growthrate2.txt"


df = pd.read_csv(dataset_txt, delimiter="\t")
df_exp = pd.read_csv(exp_dataset_txt, delimiter="\t")
df_expid = pd.read_csv(expid_dataset_txt, delimiter="\t")
df_expid_testannot = pd.read_csv(expid_dataset_testannot_txt, delimiter="\t")
df_expid_testannot_growth = pd.read_csv(expid_dataset_testannot_growth_txt, delimiter="\t")


names = set(df['Gene name'])

# drop duplicates from dataframe
df = df.drop_duplicates(['Gene name'])

# delete column 'id'
# source: https://stackoverflow.com/questions/13411544/delete-column-from-pandas-dataframe
del df['Geneid/Contrast_id']

# delete column 'LocusTag'
del df['LocusTag']
# print(df.shape)

# histogramm über Daten
# df.hist('30')
# plt.show()

#Überschriften
versuche = df_expid.iloc[0:5]



# to get what they did in their tests
# source:https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/

# new data frame with split value columns
new = df_exp["RefAnnotation"].str.split(":", n=1, expand=True)
# making separate 'Versuch' column from new data frame
df_exp["Versuch"] = new[0]
# making separate 'Anzahl' column from new data frame
df_exp["Anzahl"] = new[1]
# Dropping old Name columns
df_exp.drop(columns=["RefAnnotation"], inplace=True)
# print the data


# get only the 'Versuch' column
What_they_did_in_the_test = df_exp['Versuch']

# 'versuche' as header
versuchx = list(versuche.columns.values)
genname = df.iloc[:, 0]
#print(genname)
#print(len(genname))

# zeilen und spalten beschriftung
df.columns = versuchx
df.index = genname
del df['Genname']

# NaN to 0
df.fillna(0, inplace=True)

# alles größer 1.5 auf 1.5 setzen alles kleiner -1.5 auf -1.5 setzen
df = df.astype(float)

df[df > 1.5] = 1.5
df[df < -1.5] = -1.5

# histogramm über einen versuch
# df.hist('E-TABM-103.14.ch1-vs-E-TABM-103.3.ch1')
# plt.show()

# df.hist()
# plt.show()

'''
nach genen geclustert
# re-sort clusters by mean cluster expression levels
# nach mittelwert neu clustern

'''

# https://datatofish.com/k-means-clustering-python/

# Name gegen versuche
versuche_ohnegenname = versuchx[1:4078]
print(len(versuche_ohnegenname))


df_kmeans = pd.DataFrame(df, columns=versuche_ohnegenname)


kmeans = KMeans(n_clusters=3).fit(df_kmeans)
#centroids = kmeans.cluster_centers_


# source: https://stackoverflow.com/questions/29799053/how-to-print-result-of-clustering-in-sklearn
# labels contain all cluster numbers for rows
labels = kmeans.predict(df_kmeans)

# labels contain all cluster numbers for rows
#print(labels)

# add labels as column to dataframe
df['labels'] = labels

df_sorted_by_labels = df.sort_values('labels')


one, two, three = [], [], []
for index, row in enumerate(df_sorted_by_labels[1:-1].values):
    if df_sorted_by_labels['labels'][index] == 0:
        one.append(row[:-1])
    if df_sorted_by_labels['labels'][index] == 1:
        two.append(row[:-1])
    if df_sorted_by_labels['labels'][index] == 2:
        three.append(row[:-1])


# concat all rows to one new list (nparray)
concat_mean_one = np.concatenate(one[:], axis=0)
concat_mean_two = np.concatenate(two[:], axis=0)
concat_mean_three = np.concatenate(three[:], axis=0)


# get the means from all rows given the cluster
mean_one = np.mean(concat_mean_one)
mean_two = np.mean(concat_mean_two)
mean_three = np.mean(concat_mean_three)


one_mean = []
two_mean = []
three_mean = []
cluster_sorted_after_mean = []
for label_index, lab in enumerate(labels.tolist()):
    if df_sorted_by_labels['labels'][label_index]==0:
        if mean_one <= -0.041:
            cluster_sorted_after_mean.append(0)
            one_mean.append(0)
        if mean_one >= 0.033:
            cluster_sorted_after_mean.append(2)
            three_mean.append(2)
        if not mean_one >= 0.033 and not mean_one <= -0.041:
            cluster_sorted_after_mean.append(1)
            two_mean.append(1)
    if df_sorted_by_labels['labels'][label_index]==1:      
        if mean_two <= -0.041:
            cluster_sorted_after_mean.append(0)
            one_mean.append(0)
        if mean_two >= 0.033:
            cluster_sorted_after_mean.append(2)
            three_mean.append(2)
        if not mean_two >= 0.033 and not mean_two <= -0.041:
            cluster_sorted_after_mean.append(1)
            two_mean.append(1)
    if df_sorted_by_labels['labels'][label_index]==2:       
        if mean_three <= -0.041:
            cluster_sorted_after_mean.append(0)
            one_mean.append(0)
        if mean_three >= 0.033:
            cluster_sorted_after_mean.append(2)
            three_mean.append(2)
        if not mean_three >= 0.033 and not mean_three <= -0.041:
            cluster_sorted_after_mean.append(1)
            two_mean.append(1)


# add labels as column to dataframe
df_sorted_by_labels['labels_mean'] = cluster_sorted_after_mean

df_sorted_by_mean = df_sorted_by_labels.sort_values('labels_mean')



'''

#plt.scatter(df_kmeans['E-TABM-103.14.ch1-vs-E-TABM-103.3.ch1'], df_kmeans['E-TABM-103.19.ch1-vs-E-TABM-103.3.ch1'], c= kmeans.labels_.astype(float), s=50, alpha=0.5)
#plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
#plt.show()

'''
'''
nach experimenten geclustert
##nach größe sortieren
'''

# https://stackoverflow.com/questions/16822066/rotate-pandas-dataframe-90-degrees

# df gedreht
df_gedreht = df.T

df_kmeans_gedreht = pd.DataFrame(df_gedreht, columns=genname)


kmeans_gedreht = KMeans(n_clusters=3).fit(df_kmeans_gedreht)
#centroids_gedreht = kmeans_gedreht.cluster_centers_


# source: https://stackoverflow.com/questions/29799053/how-to-print-result-of-clustering-in-sklearn


# labels contain all cluster numbers for rows
labels_gedreht = kmeans_gedreht.predict(df_kmeans_gedreht)


# labels contain all cluster numbers for rows

# add labels as column to dataframe
df_gedreht['labels_gedreht'] = labels_gedreht


one_s, two_s, three_s = [], [], []
for index, row in enumerate(df_gedreht[1:-1].values):
    if df_gedreht['labels_gedreht'][index] == 0:
        one_s.append(row[:-1])
    if df_gedreht['labels_gedreht'][index] == 1:
        two_s.append(row[:-1])
    if df_gedreht['labels_gedreht'][index] == 2:
        three_s.append(row[:-1])


# get the means from all rows given the cluster
size_one = len(one_s)
size_two = len(two_s)
size_three = len(three_s)


cluster_sorted_after_size = []
for label_index, lab in enumerate(labels_gedreht.tolist()):
    if df_gedreht['labels_gedreht'][label_index]==0:
        if size_one <= 540:
            cluster_sorted_after_size.append(0)
        if size_one >= 1000:
            cluster_sorted_after_size.append(2)
        if not size_one >= 1000 and not size_one <= 540:
            cluster_sorted_after_size.append(1)
    if df_gedreht['labels_gedreht'][label_index]==1:      
        if size_two <= 540:
            cluster_sorted_after_size.append(0)
        if size_two >= 1000:
            cluster_sorted_after_size.append(2)
        if not size_two >= 1000 and not size_two <= 540:
            cluster_sorted_after_size.append(1)
    if df_gedreht['labels_gedreht'][label_index]==2:       
        if size_three <= 540:
            cluster_sorted_after_size.append(0)
        if size_three >= 1000:
            cluster_sorted_after_size.append(2)
        if not size_three >= 1000 and not size_three <= 540:
            cluster_sorted_after_size.append(1)


# add labels as column to dataframe
df_gedreht['labels_size'] = cluster_sorted_after_size

#https://chrisalbon.com/python/data_wrangling/pandas_dropping_column_and_rows/

#letze Spalte gelöscht da dort labels drin steht
df_gedreht_ohne_labels = df_gedreht.drop('labels')


#made labelsize to list
Versuchsgroesse = df_gedreht['labels_size']
Versuchsgroesse_Liste = Versuchsgroesse.tolist()


#löschen des letzten werts aus der liste, da sie lable enthält
Versuchsgroesse_Liste_weniger1 = Versuchsgroesse_Liste[:-1]


# eine zahl an die Versuchsliste fügen, da sonst zu wenig zahlenvorhanden sind wegen der labels zeile
Versuchsgroesse_Liste_plus1 = [Versuchsgroesse_Liste_weniger1,[50,50]]
concat_Versuchsgroesse = np.concatenate(Versuchsgroesse_Liste_plus1[:], axis=0)


#versuchsgröße als reihe an die nach mean sortierten df anhängen
df_sorted_by_mean.loc['Versuchsgroesse'] = concat_Versuchsgroesse


#den bereits nach mittelwert sortierten nach größe sortieren
df_sorted_by_mean_and_size = df_sorted_by_mean.iloc[:, np.argsort(df_sorted_by_mean.loc['Versuchsgroesse'])]


df_originalexpression = df_sorted_by_mean_and_size.apply(lambda x : 500 * 2**(x))


df_final = df_sorted_by_mean_and_size.iloc[0:4274,0:4077]
df_final_originalexpression = df_originalexpression.iloc[0:4274,0:4077]




'''
Daten für wachstumsversuche
'''

df_growthrate = df_sorted_by_mean_and_size[['GSM118600.ch1-vs-GSM118601.ch1', 'GSM118602.ch1-vs-GSM118601.ch1', 'GSM118603.ch1-vs-GSM118601.ch1', 'GSM118604.ch1-vs-GSM118601.ch1',
'GSM118605.ch1-vs-GSM118601.ch1', 'GSM118606.ch1-vs-GSM118601.ch1', 'GSM118607.ch1-vs-GSM118601.ch1', 'GSM118608.ch1-vs-GSM118601.ch1',
'GSM118609.ch1-vs-GSM118601.ch1', 'GSM118610.ch1-vs-GSM118601.ch1', 'GSM118611.ch1-vs-GSM118601.ch1', 'GSM118612.ch1-vs-GSM118601.ch1',
'GSM118613.ch1-vs-GSM118601.ch1','GSM118614.ch1-vs-GSM118601.ch1', 'GSM118615.ch1-vs-GSM118601.ch1', 'GSM118616.ch1-vs-GSM118601.ch1',
'GSM118617.ch1-vs-GSM118601.ch1', 'GSM118618.ch1-vs-GSM118601.ch1', 'GSM118619.ch1-vs-GSM118601.ch1', 'GSM118620.ch1-vs-GSM118601.ch1',
'GSM118621.ch1-vs-GSM118601.ch1', 'GSM118622.ch1-vs-GSM118601.ch1', 'GSM118623.ch1-vs-GSM118601.ch1', 'GSM702193.ch1-vs-GSM702192.ch1',
'GSM702194.ch1-vs-GSM702192.ch1', 'GSM702195.ch1-vs-GSM702192.ch1', 'GSM702196.ch1-vs-GSM702192.ch1', 'GSM702197.ch1-vs-GSM702196.ch1',
'GSM702198.ch1-vs-GSM702196.ch1', 'GSM702208.ch1-vs-GSM702196.ch1', 'GSM702209.ch1-vs-GSM702212.ch1', 'GSM702210.ch1-vs-GSM702212.ch1',
'GSM702211.ch1-vs-GSM702212.ch1', 'GSM702212.ch1-vs-GSM702192.ch1', 'GSM1248718.ch1-vs-GSM1248714.ch1', 'GSM1248715.ch1-vs-GSM1248714.ch1',
'GSM1248716.ch1-vs-GSM1248714.ch1', 'GSM1248719.ch1-vs-GSM1248714.ch1', 'GSM1248717.ch1-vs-GSM1248714.ch1', 'GSM589790.ch1-vs-GSM589790.ch2',
'GSM589791.ch1-vs-GSM589791.ch2', 'GSM589793.ch1-vs-GSM589793.ch2', 'GSM589792.ch1-vs-GSM589792.ch2', 'GSM589794.ch1-vs-GSM589794.ch2',
'GSM589796.ch1-vs-GSM589796.ch2', 'GSM589795.ch1-vs-GSM589795.ch2','labels_mean']]


#growthrate nach groesse sortieren
df_growthrate_sorted_by_size = df_growthrate.iloc[:, np.argsort(df_growthrate.loc['Versuchsgroesse'])]


df_growthrate_originalexpression = df_growthrate_sorted_by_size.apply(lambda x : 500 * 2**(x))


#für heatmap anpassten salso label und größe raus schmeißen
df_growthrate_final = df_growthrate_sorted_by_size.iloc[0:4274,0:46]


df_growthrate_originalespression_final = df_growthrate_originalexpression.iloc[0:4274,0:46]


#plt.scatter(df_gedreht['thrL'], df_gedreht['thrA'], c= kmeans.labels_.astype(float), s=50, alpha=0.5)
#plt.scatter(centroids_gedreht[:, 0], centroids_gedreht[:, 1], c='red', s=50)
# plt.show()



'''
sortieren nach expressionsleveln 
'''

#gibt einen dataframe mit nur dem ersten cluster aus
df_erstes_cluster = df_sorted_by_mean_and_size.iloc[0:len(one),:]
df_zweites_cluster = df_sorted_by_mean_and_size.iloc[len(one):len(one)+len(two),:]
df_drittes_cluster = df_sorted_by_mean_and_size.iloc[len(one)+len(two):len(one)+len(two)+len(three),:]

print(len(df_erstes_cluster))
print(len(df_zweites_cluster))
print(len(df_drittes_cluster))
print(len(df_erstes_cluster)+len(df_drittes_cluster)+len(df_zweites_cluster))

#means aller Versuche sortiert nach cluster 1 und dann in  liste gespeichert
df_mean_cluster1_Versuche = df_erstes_cluster.mean(axis = 0)
df_mean_cluster1_list = df_mean_cluster1_Versuche.tolist()


#mean von cluster1 als Reihe an df der nach mean und Groeße sortiert ist gehängt
df_sorted_by_mean_and_size.loc['mean_cluster1'] = df_mean_cluster1_list


# nach der zuvor erstellten Reihe Sortieren
df_sorted_by_meancluster1 = df_sorted_by_mean_and_size.iloc[:, np.argsort(df_sorted_by_mean_and_size.loc['mean_cluster1'])]


#df nach groesse der cluster sortieren
df_sorted_by_size = df_gedreht.sort_values('labels_size')


#gibt einen dataframe mit nur dem ersten cluster aus
df_erstes_cluster_gene = df_sorted_by_size.iloc[0:len(one_s),:]


#means aller gene sortieren nach cluster 1 und dann als Liste gespeichert
df_mean_cluster1_gene = df_erstes_cluster_gene.mean(axis = 0)
df_mean_cluster1_gene_list = df_mean_cluster1_gene.tolist()


#mean von cluster1gene als Spalte an df der nach meancluster1 und Groeße sortiert ist gehängt
df_sorted_by_meancluster1['mean_cluster1_gene'] = df_mean_cluster1_gene_list


#deleten unerwünschter Zeilen und Spalten
df_cluster1_0 = df_sorted_by_meancluster1.drop(['Versuchsgroesse', 'mean_cluster1'])

df_cluster1_1 = df_cluster1_0.drop('labels', axis=1)
df_cluster1_2 = df_cluster1_1.drop('labels_mean', axis=1)
df_cluster1_final = df_cluster1_2.drop('mean_cluster1_gene', axis=1)




'''
nur Wachstumsversuche sortieren nach expression
'''
df_growthrate_expression = df_sorted_by_meancluster1[['GSM118600.ch1-vs-GSM118601.ch1', 'GSM118602.ch1-vs-GSM118601.ch1', 'GSM118603.ch1-vs-GSM118601.ch1', 'GSM118604.ch1-vs-GSM118601.ch1',
'GSM118605.ch1-vs-GSM118601.ch1', 'GSM118606.ch1-vs-GSM118601.ch1', 'GSM118607.ch1-vs-GSM118601.ch1', 'GSM118608.ch1-vs-GSM118601.ch1',
'GSM118609.ch1-vs-GSM118601.ch1', 'GSM118610.ch1-vs-GSM118601.ch1', 'GSM118611.ch1-vs-GSM118601.ch1', 'GSM118612.ch1-vs-GSM118601.ch1',
'GSM118613.ch1-vs-GSM118601.ch1','GSM118614.ch1-vs-GSM118601.ch1', 'GSM118615.ch1-vs-GSM118601.ch1', 'GSM118616.ch1-vs-GSM118601.ch1',
'GSM118617.ch1-vs-GSM118601.ch1', 'GSM118618.ch1-vs-GSM118601.ch1', 'GSM118619.ch1-vs-GSM118601.ch1', 'GSM118620.ch1-vs-GSM118601.ch1',
'GSM118621.ch1-vs-GSM118601.ch1', 'GSM118622.ch1-vs-GSM118601.ch1', 'GSM118623.ch1-vs-GSM118601.ch1', 'GSM702193.ch1-vs-GSM702192.ch1',
'GSM702194.ch1-vs-GSM702192.ch1', 'GSM702195.ch1-vs-GSM702192.ch1', 'GSM702196.ch1-vs-GSM702192.ch1', 'GSM702197.ch1-vs-GSM702196.ch1',
'GSM702198.ch1-vs-GSM702196.ch1', 'GSM702208.ch1-vs-GSM702196.ch1', 'GSM702209.ch1-vs-GSM702212.ch1', 'GSM702210.ch1-vs-GSM702212.ch1',
'GSM702211.ch1-vs-GSM702212.ch1', 'GSM702212.ch1-vs-GSM702192.ch1', 'GSM1248718.ch1-vs-GSM1248714.ch1', 'GSM1248715.ch1-vs-GSM1248714.ch1',
'GSM1248716.ch1-vs-GSM1248714.ch1', 'GSM1248719.ch1-vs-GSM1248714.ch1', 'GSM1248717.ch1-vs-GSM1248714.ch1', 'GSM589790.ch1-vs-GSM589790.ch2',
'GSM589791.ch1-vs-GSM589791.ch2', 'GSM589793.ch1-vs-GSM589793.ch2', 'GSM589792.ch1-vs-GSM589792.ch2', 'GSM589794.ch1-vs-GSM589794.ch2',
'GSM589796.ch1-vs-GSM589796.ch2', 'GSM589795.ch1-vs-GSM589795.ch2','mean_cluster1_gene']]


# nach der zuvor erstellten Reihe Sortieren
df_sorted_meancluster1_expression = df_growthrate_expression.iloc[:, np.argsort(df_growthrate_expression.loc['mean_cluster1'])]


#für heatmap anpassten salso label und größe raus schmeißen
df_growthrate_expression_final = df_sorted_meancluster1_expression.iloc[0:4274,0:46]


'''
Liniendiagramm 
'''

#versuchsüberschriften als reihenbeschriftung zur besseren übersicht
Versuche_doppelt = df_expid_testannot.iloc[:, 0]
df_expid_testannot.index = Versuche_doppelt


#nur die versuche für Growthrate in df speichern
neue_df = df_expid_testannot.loc[['GSM118600.ch1-vs-GSM118601.ch1', 'GSM118602.ch1-vs-GSM118601.ch1', 'GSM118603.ch1-vs-GSM118601.ch1', 'GSM118604.ch1-vs-GSM118601.ch1',
'GSM118605.ch1-vs-GSM118601.ch1', 'GSM118606.ch1-vs-GSM118601.ch1', 'GSM118607.ch1-vs-GSM118601.ch1', 'GSM118608.ch1-vs-GSM118601.ch1',
'GSM118609.ch1-vs-GSM118601.ch1', 'GSM118610.ch1-vs-GSM118601.ch1', 'GSM118611.ch1-vs-GSM118601.ch1', 'GSM118612.ch1-vs-GSM118601.ch1',
'GSM118613.ch1-vs-GSM118601.ch1','GSM118614.ch1-vs-GSM118601.ch1', 'GSM118615.ch1-vs-GSM118601.ch1', 'GSM118616.ch1-vs-GSM118601.ch1',
'GSM118617.ch1-vs-GSM118601.ch1', 'GSM118618.ch1-vs-GSM118601.ch1', 'GSM118619.ch1-vs-GSM118601.ch1', 'GSM118620.ch1-vs-GSM118601.ch1',
'GSM118621.ch1-vs-GSM118601.ch1', 'GSM118622.ch1-vs-GSM118601.ch1', 'GSM118623.ch1-vs-GSM118601.ch1', 'GSM702193.ch1-vs-GSM702192.ch1',
'GSM702194.ch1-vs-GSM702192.ch1', 'GSM702195.ch1-vs-GSM702192.ch1', 'GSM702196.ch1-vs-GSM702192.ch1', 'GSM702197.ch1-vs-GSM702196.ch1',
'GSM702198.ch1-vs-GSM702196.ch1', 'GSM702208.ch1-vs-GSM702196.ch1', 'GSM702209.ch1-vs-GSM702212.ch1', 'GSM702210.ch1-vs-GSM702212.ch1',
'GSM702211.ch1-vs-GSM702212.ch1', 'GSM702212.ch1-vs-GSM702192.ch1', 'GSM1248718.ch1-vs-GSM1248714.ch1', 'GSM1248715.ch1-vs-GSM1248714.ch1',
'GSM1248716.ch1-vs-GSM1248714.ch1', 'GSM1248719.ch1-vs-GSM1248714.ch1', 'GSM1248717.ch1-vs-GSM1248714.ch1', 'GSM589790.ch1-vs-GSM589790.ch2',
'GSM589791.ch1-vs-GSM589791.ch2', 'GSM589793.ch1-vs-GSM589793.ch2', 'GSM589792.ch1-vs-GSM589792.ch2', 'GSM589794.ch1-vs-GSM589794.ch2',
'GSM589796.ch1-vs-GSM589796.ch2', 'GSM589795.ch1-vs-GSM589795.ch2']]


###########neue_df.to_csv(r'growthrate.txt', header=None, index=None, sep=' ', mode='a')


# überschriften der Spalten für grothrate df
df_expid_testannot_growth.columns = ['Versuche','Growth_Rate']


#row aus df als liste speichern um an df_expid_testannot_growth als spalte anzuhängen
meancluster1_row = df_growthrate_expression.iloc[-1]
meancluster1_liste = meancluster1_row.tolist()
meancluster1_ohneletzterwert = meancluster1_liste[0:46]


#meancluster1_liste als Spalte an df_expid_testannot_growth  angehängt
df_expid_testannot_growth['meancluster1_versuche'] = meancluster1_ohneletzterwert


#row aus df als liste speichern um an df_expid_testannot_growth als spalte anzuhängen (zur sortierung vach versuchsgröße) 
Versuchsgroesse_row = df_growthrate.iloc[-1]
Versuchsgroesse_liste = Versuchsgroesse_row.tolist()
Versuchsgroesse_ohneletzterwert = Versuchsgroesse_liste[0:46]


#meancluster1_liste als Spalte an df_expid_testannot_growth  angehängt
df_expid_testannot_growth['Versuchsgroesse'] = Versuchsgroesse_ohneletzterwert



df_growthrate_sorted_Versuchsgroesse = df_expid_testannot_growth.sort_values('Versuchsgroesse')


#in cluster gesplitet
df1 = df_growthrate_sorted_Versuchsgroesse.iloc[0:10]
df2 = df_growthrate_sorted_Versuchsgroesse.iloc[10:45]


# überschriften der Spalten für grothrate df nach Growthrate sortiert
df1.columns = ['Versuche','Growth_Rate','meancluster1_versuche', 'Versuchsgroesse']


#gesplitete cluster sortieren nach Zeit
df1_sorted_time = df1.sort_values('Growth_Rate')
df2_sorted_time = df2.sort_values('Growth_Rate')

'''
plt.plot( 'Versuche', 'Versuchsgroesse', data=df1_sorted_time, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
plt.plot( 'Versuche', 'Versuchsgroesse', data=df2_sorted_time, marker='', color='olive', linewidth=2)
plt.legend()
plt.show()
'''

# überschriften der Spalten für grothrate df nach Growthrate sortiert
df_expid_testannot.columns = ['ContrastName','TestAnnotation']


df_expo1 = df_expid_testannot[df_expid_testannot.TestAnnotation == 'GROWTH.EXPONENTIAL:1']
df_expoplus1 = df_expid_testannot[df_expid_testannot.TestAnnotation == 'GROWTH.EXPONENTIAL:+1']
df_stat1 = df_expid_testannot[df_expid_testannot.TestAnnotation == 'GROWTH.STATIONARY:1']
df_statplus1 = df_expid_testannot[df_expid_testannot.TestAnnotation == 'GROWTH.STATIONARY:+1']



listAnnotation = df_expid_testannot.TestAnnotation.tolist()


zustaende = []
ContrastName_list = []
for label_index, lab in enumerate(listAnnotation):
    if listAnnotation[label_index] == 'GROWTH.EXPONENTIAL:1':
        zustaende.append("1")
        ContrastName_list.append(df_expid_testannot["ContrastName"][label_index])
    if listAnnotation[label_index] == 'GROWTH.EXPONENTIAL:+1':
        zustaende.append("1")
        ContrastName_list.append(df_expid_testannot["ContrastName"][label_index])
    if listAnnotation[label_index] == 'GROWTH.STATIONARY:1':
        zustaende.append("2")
        ContrastName_list.append(df_expid_testannot["ContrastName"][label_index])
    if listAnnotation[label_index] == 'GROWTH.STATIONARY:+1':
        zustaende.append("2")
        ContrastName_list.append(df_expid_testannot["ContrastName"][label_index])
        


#beide listen zu einem DataFrame zusammenführen
d = {'Versuche':ContrastName_list,'Growth_condition':zustaende}
df_growth_contions = pd.DataFrame(d, columns=['Versuche','Growth_condition'])


#versuche zur Liste
Versuche_growth = df_growth_contions["Versuche"].tolist()


#nur die versuche für Growth in df speichern
growth_df = df_gedreht.loc[Versuche_growth]

# liste der Versuche mit condition
list_Versuche_cluster = growth_df['labels_size'].values

#lister der cluster an die condition liste hinzufügen
df_growth_contions['Cluster'] = list_Versuche_cluster


# sortieren nach cluster column
df_growth_contions_sorted_cluster = df_growth_contions.sort_values('Cluster')


# Versuche in cluster einteilen 
ones_cluster = []
two_cluster = []
three_cluster = []

for label_index, lab in enumerate(list_Versuche_cluster):
    if df_growth_contions_sorted_cluster["Cluster"][label_index] == 0:
        ones_cluster.append("1")
    if df_growth_contions_sorted_cluster["Cluster"][label_index] == 1:
        two_cluster.append("2")
    if df_growth_contions_sorted_cluster["Cluster"][label_index] == 2:
        three_cluster.append("3")
    

df_erstes_cluster_growth = df_growth_contions_sorted_cluster.iloc[0:(len(ones_cluster)),:]
df_zweite_cluster_growth = df_growth_contions_sorted_cluster.iloc[len(ones_cluster):(len(ones_cluster)+len(two_cluster)),:]
df_drittes_cluster_growth = df_growth_contions_sorted_cluster.iloc[(len(ones_cluster)+len(two_cluster)):(len(ones_cluster)+len(two_cluster)+len(three_cluster)),:]



'''
# multiple line plot
plt.plot( 'Versuche', 'Growth_condition', data=df_erstes_cluster_growth, marker='o', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4)
plt.plot( 'Versuche', 'Growth_condition', data=df_zweite_cluster_growth, marker='', color='olive', linewidth=2)
plt.plot( 'Versuche', 'Growth_condition', data=df_drittes_cluster_growth, marker='', color='olive', linewidth=2, linestyle='dashed', label="toto")
plt.legend()
plt.show()

'''



# liste der Gene mit 
list_Gene_cluster1 = df_cluster1_0['labels'].values


# Gene in cluster einteilen 
ones_cluster_gene = []
two_cluster_gene = []
three_cluster_gene = []

for label_index, lab in enumerate(list_Gene_cluster1):
    if df_cluster1_0["labels_mean"][label_index] == 0:
        ones_cluster_gene.append("1")
    if df_cluster1_0["labels_mean"][label_index] == 1:
        two_cluster_gene.append("2")
    if df_cluster1_0["labels_mean"][label_index] == 2:
        three_cluster_gene.append("3")



#ein df zu dreien nach cluster geteilt
df_erstes_cluster1 = df_cluster1_0.iloc[0:(len(ones_cluster_gene)),:]
df_zweite_cluster1 = df_cluster1_0.iloc[len(ones_cluster_gene):(len(ones_cluster_gene)+len(two_cluster_gene)),:]
df_drittes_cluster1 = df_cluster1_0.iloc[(len(ones_cluster_gene)+len(two_cluster_gene)):(len(ones_cluster_gene)+len(two_cluster_gene)+len(three_cluster_gene)),:]


gene_in_cluster_1 = list(df_erstes_cluster1.index.values)
gene_in_cluster_2 = list(df_zweite_cluster1.index.values)
gene_in_cluster_3 = list(df_drittes_cluster1.index.values)

#print(len(gene_in_cluster_1)+len(gene_in_cluster_2)+len(gene_in_cluster_3))
#print(len(gene_in_cluster_1))
#print(len(gene_in_cluster_2))
#print(len(gene_in_cluster_3))

#summe für expression jedes clusters einzeln
expression_cluster1 = df_erstes_cluster1.mean(axis=0) 
expression_cluster2 = df_zweite_cluster1.mean(axis=0) 
expression_cluster3 = df_drittes_cluster1.mean(axis=0) 

#summe für expression jedes clusters einzeln als Liste
expression_cluster1_list = expression_cluster1.tolist()
expression_cluster2_list = expression_cluster2.tolist()
expression_cluster3_list = expression_cluster3.tolist()


#Cluster als Reihe an den df_cluster1_0 
df_cluster1_0.loc['cluster1'] = expression_cluster1_list
df_cluster1_0.loc['cluster2'] = expression_cluster2_list
df_cluster1_0.loc['cluster3'] = expression_cluster3_list


df_cluster1_1_neu = df_cluster1_0.drop('labels', axis=1)
df_cluster1_2_neu = df_cluster1_1_neu.drop('labels_mean', axis=1)
df_cluster1_final_neu = df_cluster1_2_neu.drop('mean_cluster1_gene', axis=1)


#cluster drehen
df_cluster1_gedreht = df_cluster1_final_neu.T


#cluster ueberschriften
columnueberschrift = list(df_cluster1_final_neu.index.values)
df_cluster1_gedreht.columns = columnueberschrift

indexueberschrift = list(df_cluster1_final_neu.columns.values)


new = df_cluster1_gedreht[['cluster1', 'cluster2', 'cluster3']].copy()

# add Versuche as column to dataframe
new['Versuche'] = indexueberschrift

'''
# multiple line plot
plt.figure(figsize=(10,10))
plt.plot( 'Versuche', 'cluster3', data=new, color='skyblue', linewidth=2)
#plt.plot( 'Versuche', 'cluster2', data=new, color='olive', linewidth=2)
plt.plot( 'Versuche', 'cluster1', data=new, color='red', linewidth=2)
plt.xticks(rotation=90)
plt.title("durchschnittliche Clusterexpression für die sortierten Experimente")
plt.xlabel("Versuche")
plt.ylabel("durchschnittliche Expression")
plt.legend()
plt.show()
'''
'''
# multiple line plot
plt.figure(figsize=(20,10))
plt.plot( 'Versuche', 'cluster3', data=new, color='olive', linewidth=2,label='Transport- und metabolische Enzyme-Cluster')
plt.plot( 'Versuche', 'cluster2', data=new, color='skyblue', linewidth=2, label='Wachstumsunabhängige Proteine-Cluster')
plt.plot( 'Versuche', 'cluster1', data=new, color='red', linewidth=2, label='Ribosomen-Cluster')
test = ["Transport- und metabolische Enzyme", "Wachstumsunabhängige Proteine", "Ribosomen"]
plt.legend(loc='lower center',prop={'size': 15})
#plt.title("durchschnittliche Clusterexpression für die sortierten Experimente")
plt.xlabel("Experimente", size = 20)
plt.ylabel("durchschnittliche Expression", size = 20)
plt.xticks([])
plt.show()
'''



# multiple line plot
plt.figure(figsize=(10,10))
plt.scatter( 'cluster1', 'cluster3', data=new, color='skyblue')
#plt.plot( 'Versuche', 'cluster2', data=new, color='olive', linewidth=2)
plt.xlabel("ribosomale Proteine", size = 15)
plt.ylabel("metabolische und Transportenzyme", size = 15)
plt.show()



'''
##########
versuche zur korreltaion
##########

# Get dataset
pr = [24.962410217100825, 24.68808625174456, 24.413504712744786, 24.138736997770838, 23.863856707564114, 23.58893982182796, 23.314064902511983, 23.03931332918861, 22.764769572102665, 22.490521509621775, 22.21666079816729, 21.943283304423314, 21.67048961173044, 21.39838561521231, 21.127083223515463, 20.856701189224726, 20.587366095385317, 20.319213532368277, 20.052389508164183, 19.78705214663441, 19.523373743241855, 19.261543267584937, 19.001769428422705, 18.744284452324703, 18.489348775204576, 18.23725691211455, 17.98834486245842, 17.742999536854473, 17.501670875931932, 17.264887597050645, 17.03327789353109, 16.807596980967467, 16.588764193713747, 16.3779132597351, 16.176458963090756, 15.986171535326488, 15.809161591450145, 15.647087075617883, 15.494907024033026]
pt = [7398.976343780549, 7446.309170515266, 7494.395058224955, 7543.259006693594, 7592.9275845885195, 7643.429074749881, 7694.793636385929, 7747.0534865173095, 7800.243103392513, 7854.399455050469, 7909.562256745338, 7965.774261595364, 8023.081589596188, 8081.534101079319, 8141.185821837949, 8202.095428534567, 8264.32680471376, 8327.949679842906, 8393.040366419573, 8459.6826134292, 8527.968598523008, 8598.000086438891, 8669.88978774964, 8743.762960428434, 8819.759307587805, 8898.035238910177, 8978.766581915568, 9062.151853940475, 9148.416238892089, 9237.816457687439, 9330.6467819293, 9427.246519229822, 9528.009378739232, 9633.39505978259, 9743.94224325276, 9860.273847034547, 9983.030492523558, 10112.302504634892, 10243.584865032532]


def slope_intercept(x_val,y_val):
    x = np.array(x_val)
    y = np.array(y_val)
    m = (((np.mean(x)*np.mean(y)) - np.mean(x*y))/((np.mean(x)*np.mean(x)) - np.mean(x*x)))
    m = round(m,2)
    b = (np.mean(y) - np.mean(x)*m)
    b = round(b,2)

    return m,b
print(slope_intercept(pr,pt))

m,b = slope_intercept(pr,pt)

reg_line = [(m*x)+b for x in pr]
print(reg_line)
print(len(reg_line))


Raten_df = pd.DataFrame(
    {'pr': pr,
     'pt': pt,
     'rl': reg_line
    })




# Normalize total_bedrooms column
x_array = np.array(Raten_df['pr'])
y_array = np.array(Raten_df['pt'])
rl_array = np.array(Raten_df['rl'])
cl1_array = np.array(new['cluster1'])
cl2_array = np.array(new['cluster3'])
normalized_X = preprocessing.normalize([x_array])
normalized_Y = preprocessing.normalize([y_array])
normalized_RL = preprocessing.normalize([rl_array])
normalized_cl1 = preprocessing.normalize([cl1_array])
normalized_cl2 = preprocessing.normalize([cl2_array])



cluster1_list = new["cluster1"].tolist()
cluster2_list = new["cluster2"].tolist()
#print(cluster1_list)



myIntpr = np.max(pr)
newListpr = [x / myIntpr for x in pr]
#print(newListpr)

myIntpt = np.max(pt)
newListpt = [x / myIntpt for x in pt]
#print(newListpt)

myIntcl1 = np.max(cluster1_list)
newListcl1 = [x / myIntcl1 for x in cluster1_list]
#print(newListcl1)

myIntcl2 = np.max(cluster2_list)
newListcl2 = [x / myIntcl2 for x in cluster2_list]
#print(newListcl2)



#import scipy

#slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(pr, pt)
#print(slope)
#print(intercept)
#print(r_value)
#print(p_value)



myIntrl = np.max(reg_line)
newListrl = [x / myIntrl for x in reg_line]
'''

'''
plt.scatter(pr,pt, color='red')
plt.plot(pr,reg_line)
plt.ylabel('bla')
plt.xlabel('blabla')
plt.title('lineare regression')
plt.show()
'''
'''
# multiple line plot
plt.figure(figsize=(10,10))
plt.scatter( 'cluster1', 'cluster3', data=new, color='skyblue')
plt.scatter( pr, pt, color='red')
#plt.scatter( newListcl1, newListcl2, color='green')
plt.plot(normalized_RL, color='green')
plt.xlabel("cluster1")
plt.ylabel("cluster3")
plt.show()

'''



'''
# multiple line plot
plt.figure(figsize=(10,10))
plt.scatter( 'cluster1', 'cluster3', data=new, color='skyblue')
#plt.scatter( normalized_X, normalized_Y, color='red')
#plt.scatter( normalized_cl1, normalized_cl2, color='green')
plt.xlabel("ribosomale Proteine")
plt.ylabel("metabolische und Transportenzyme")
plt.show()
'''



df_sorted_my_mean_and_size_gedreht = df_sorted_by_mean_and_size.T


geeene = df_final.index.values
geeene = geeene.tolist()
#print(geeene)

#print(df_sorted_by_mean_and_size)
#print(df_sorted_by_mean_and_size.iloc[0:-2])

# Versuche in cluster einteilen 
ones_gene = []
two_gene = []
three_gene = []

for label_index, lab in enumerate(df_sorted_by_mean_and_size.iloc[0:-2]):
    if df_sorted_by_mean_and_size["labels_mean"][label_index] == 0:
        ones_gene.append("1")
    if df_sorted_by_mean_and_size["labels_mean"][label_index] == 1:
        two_gene.append("2")
    if df_sorted_by_mean_and_size["labels_mean"][label_index] == 2:
        three_gene.append("3")
    
#print(len(ones_gene))
#print(len(two_gene))
#print(len(three_gene))

gene_cluster1 = geeene[0:len(ones_gene)]
gene_cluster2 = geeene[len(ones_gene):len(ones_gene)+len(two_gene)]
gene_cluster3 = geeene[len(ones_gene)+len(two_gene):len(ones_gene)+len(two_gene)+len(three_gene)]

#print(gene_cluster1)
#print(gene_cluster2)
#print(len(gene_cluster3))
        



'''
one_groesse, two_groesse, three_groesse= [], [], []
for index, row in enumerate(df_sorted_my_mean_and_size_gedreht): #[1:-1].values
    if df_sorted_by_mean_and_size['Versuchsgroesse'][index] == 0:
        one_groesse.append(row[:-1])
    if df_sorted_by_mean_and_size['Versuchsgroesse'][index] == 1:
        two_groesse.append(row[:-1])
    if df_sorted_by_mean_and_size['Versuchsgroesse'][index] == 2:
        three_groesse.append(row[:-1])

'''

'''

#expression aller Cluster
# plot using a color palette
sns.heatmap(df_final, cmap='coolwarm', xticklabels=0, yticklabels=0)#, cmap="RdYlGn")das ist das wichtige => RdBu
#plt.title("Test") 
plt.legend(loc = 'center right',prop={'size': 30}, ncol=2, bbox_to_anchor=[0.7,1])
plt.xlabel("Experiments")
plt.ylabel("Genes")
 
#add this after your favorite color to show the plot
plt.show()
'''

'''
# plot using a color palette
sns.heatmap(df_cluster1_final, cmap='coolwarm', xticklabels=0, yticklabels=0)#, cmap="RdYlGn")das ist das wichtige => RdBu
#plt.title("Test")
plt.xlabel("Experiments")
plt.ylabel("Genes")
 
#add this after your favorite color to show the plot
plt.show()
'''