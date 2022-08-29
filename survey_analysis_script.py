# %%
import pandas as pd
df = pd.read_excel("C:\\Users\\Victoria\\Documents\\Postgraduate - University of Manchester - Genomic Medicine MSc\\Research Project\\Dissertation Final\\Variant_Classification_Survey.xlsx")

# %% [markdown]
# 

# %%
df.head()

# %%
df.iloc[:,5].value_counts()

# %%
df.loc[df.iloc[:,5] == "Bioinformatician and researcher and lecturer", "Are you completing this as a scientist, clinician, academic or other?"] = "Academic/ Researcher"

# %%
df.iloc[:,5].value_counts()

# %%
df.iloc[:,13]= df.iloc[:,13].apply(lambda x : str(x)[:-1])

# %%
df.rename(columns={cname:    str(cname.encode()).replace("\\xc2"," ").replace("\\xa0"," ")[1:].strip("'") for  cname in df.columns},inplace=True)
for c in df.columns:
    print(c.encode("utf8"))

#df

# %%
variantbatchnumber='How many variants do you normally classify in one batch?'

# %%
def sanitise_variant_batch_number (current_data):
    # if the data is not a string do something

    convert_to_actual = {
        # key : value
        "02": "Two to Five",
        "06": "Six to Ten"}

    if type(current_data) != str:
        key = current_data.strftime("%m")
        # get the key otherwise return the current_data
        return convert_to_actual.get(key,current_data)
    # always return something back otherwise you will have an empty cell in the table
    return current_data

# %%
df[variantbatchnumber] = df[variantbatchnumber].apply(sanitise_variant_batch_number)
df.head()

# %%
df = df.fillna(value=0)

# %%
pip install matplotlib

# %%
import matplotlib.pyplot as plt

# %%
genome_build = 'What genome build are you using, hg19 or hg38? Please state the patch number under other if known.'

# %%
job_title = 'Are you completing this as a scientist, clinician, academic or other?'

# %%
df[genome_build].value_counts()

# %%
df[job_title].value_counts()

# %%
gene_level_metrics = 'What gene-level metrics do you use to facilitate classification of sequence variants?'

# %%
metrics = df[gene_level_metrics].str.contains('pLI')
metrics.value_counts()

# %%
metrics = df[gene_level_metrics].str.split(';', expand=True)
metrics.head()

# %%
df[metrics].stack()

# %%
missense_variants = 'What in-silico predictors do you use to help classify missense variants?'

# %%
import textwrap
import seaborn as sns


plt.figure(figsize = (10,8))
sns.countplot(y = df.iloc[:,5],).set_yticklabels([textwrap.fill(e,16) for e in df.iloc[:,5].unique()])# for horizontal&vertical , y&x
#plt.xticks(rotation=90)  # for rotating x label 
plt.title("Profession")  # for the title
plt.ylabel('Survey Respondents Job Title')
plt.xlabel('Number of survey respondents')
plt.show()


# %%
df.loc[(df.iloc[:,7].astype(str).str.contains("ACMG")),"Guideline version"] = "ACMG"
df.loc[(df.iloc[:,7].astype(str).str.contains("ACGS")),"Guideline version"] = "ACGS"

# %%
import seaborn as sns


plt.figure(figsize = (10,5))
sns.countplot(x = df.iloc[:,-1], palette="PuBuGn")
plt.ylabel('Number of survey respondents')

plt.show()

# %%
import seaborn as sns


plt.figure(figsize = (10,15))
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,9], palette = 'PuBuGn')
plt.ylabel('Survey Respondents Job Title')
#plt.xticks(rotation=90)
plt.xlabel('Number of Respondents')
plt.show()

# %%
software_types = ["Congenica", "Varsome", "Alamut", "Golden Helix", "Variant Effect Predictor", "Exomiser", "Other (please specify below in question 6)"]

# %%
data = df.copy()

# %%

plt.figure(figsize = (20,10))
data = pd.DataFrame(data[software_types])
df1 = data.melt(var_name='Software Usage', value_name='Software Name')

sns.countplot(y='Software Usage', hue='Software Name', data=df1, palette="Paired")
plt.legend()
plt.show()

# %%
df1 = data.melt(var_name='Software Usage', value_name='Software Name')
df_plot = df1.groupby(['Software Usage', 'Software Name']).size().reset_index().pivot(columns = 'Software Name', index = 'Software Usage', values = 0)

# %%
df1

# %%
df1 = data.melt(var_name='Software Usage', value_name='Software Name')
df_plot = df1.groupby(['Software Usage', 'Software Name']).size().reset_index().pivot(columns = 'Software Name', index = 'Software Usage', values = 0)
plt.figure(figsize = (15,60))
ax = df_plot.plot(kind = 'bar', stacked = True)



# %%
from turtle import color
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot.plot(kind = 'barh', stacked = True, cmap = "YlOrRd")
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')
plt.show()

# %%
df.loc[df.Congenica.isin(["0","1","2","3","4","5"]), "congenica_user"] = 1
df.loc[df.Varsome.isin(["0","1","2","3","4","5"]), "Varsome_user"] = 1
df.loc[df.Alamut.isin(["0","1","2","3","4","5"]), "Alamut_user"] = 1
df.loc[df["Golden Helix"].isin(["0","1","2","3","4","5"]), "golden_helix_user"] = 1
df.loc[df["Variant Effect Predictor"].isin(["0","1","2","3","4","5"]), "Variant_effect_predictor_user"] = 1
df.loc[df.Exomiser.isin(["0","1","2","3","4","5"]), "Exomiser_user"] = 1
df.loc[df['Other (please specify below in question 6)'].isin(["0","1","2","3","4","5"]), "other_user"] = 1

# %%
user_df = pd.concat([df.iloc[:,5],df.iloc[:,-6:]], 1)

# %%
user_df.fillna(0, inplace = True)

# %%
user_df.head()

# %%
df1 = user_df.melt(id_vars='Are you completing this as a scientist, clinician, academic or other?')
df1 = df1.drop(df1.loc[df1.value ==0].index)
df_plot = df1.groupby(['Are you completing this as a scientist, clinician, academic or other?', 'variable']).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?', index = 'variable', values = 0)

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot.plot(kind = 'barh', stacked = True, cmap = "Spectral_r")
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    


plt.show()

# %%
import seaborn as sns
fig, ax = plt.subplots(2,3)
#sns.set(rc={'figure.figsize':(45.7,15.27)})
fig.set_figheight(15)
fig.set_figwidth(25)
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,-6], ax= ax[0][0], palette="RdPu")
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,-5], ax= ax[0][1], palette="Dark2_r")
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,-4], ax= ax[0][2], palette="cividis_r")
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,-3], ax= ax[1][0], palette="gist_stern")
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,-2], ax= ax[1][1], palette="inferno_r")
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,-1], ax= ax[1][2], palette="nipy_spectral")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# %%
import seaborn as sns


plt.figure(figsize = (10,15))
sns.countplot(y = df.iloc[:,5], hue = df.iloc[:,8], palette="PuBuGn")
plt.ylabel('Survey Respondents Job Title')
#plt.xticks(rotation=90)
plt.show()

# %%
predictors_df = df.loc[:,df.columns.astype(str).str.contains("predictors")]

# %%
predictors_df = pd.concat([df.iloc[:,5], predictors_df], 1)

# %%
import warnings
warnings.filterwarnings('ignore')

# %%
temp_df = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,1].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name = predictors_df.columns[1] )
temp_df2 = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,2].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name = predictors_df.columns[2])
temp_df3 = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,3].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name =  predictors_df.columns[3])
temp_df4 = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,4].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name = predictors_df.columns[4] )
temp_df5 = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,5].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name =  predictors_df.columns[5])
temp_df6 = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,6].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name =  predictors_df.columns[6])
temp_df7 = pd.concat([predictors_df.iloc[:,0],predictors_df.iloc[:,7].str.split(";", 7, expand = True)],1).melt(id_vars='Are you completing this as a scientist, clinician, academic or other?', var_name =  predictors_df.columns[7])

# %%
temp_df = temp_df.loc[(temp_df.value != "") & (temp_df.value != "None")]
temp_df2 = temp_df2.loc[(temp_df2.value != "") & (temp_df2.value != "None")]
temp_df3 = temp_df3.loc[(temp_df3.value != "") & (temp_df3.value != "None")]
temp_df4 = temp_df4.loc[(temp_df4.value != "") & (temp_df4.value != "None") & (temp_df4.value != 'Would like to know more about good options for this')]
temp_df5 = temp_df5.loc[(temp_df5.value != "") & (temp_df5.value != "None")]
temp_df6 = temp_df6.loc[(temp_df6.value != "") & (temp_df6.value != "None")]
temp_df7 = temp_df7.loc[(temp_df7.value != "") & (temp_df7.value != "None")]

# %%
temp_df2.loc[temp_df2.value == "Align GVD;", "value"] = "Align GVD"
temp_df3.loc[temp_df3.value == "Align GVD;", "value"] = "Align GVD"
temp_df4.loc[temp_df4.value == "Align GVD;", "value"] = "Align GVD"
temp_df5.loc[temp_df5.value == "Align GVD;", "value"] = "Align GVD"
temp_df6.loc[temp_df6.value == "Align GVD;", "value"] = "Align GVD"
temp_df7.loc[temp_df7.value == "Align GVD;", "value"] = "Align GVD"

# %%
temp_df.drop(temp_df.loc[temp_df.value == "ACMG"].index, inplace = True)
temp_df2.drop(temp_df2.loc[temp_df2.value == "ACMG"].index, inplace = True)
temp_df3.drop(temp_df3.loc[temp_df3.value == "ACMG"].index, inplace = True)
temp_df4.drop(temp_df4.loc[temp_df4.value == "ACMG"].index, inplace = True)
temp_df5.drop(temp_df5.loc[temp_df5.value == "ACMG"].index, inplace = True)
temp_df6.drop(temp_df6.loc[temp_df6.value == "ACMG"].index, inplace = True)
temp_df7.drop(temp_df7.loc[temp_df7.value == "ACMG"].index, inplace = True)

# %%
temp_df.loc[temp_df.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"
temp_df2.loc[temp_df2.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"
temp_df3.loc[temp_df3.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"
temp_df4.loc[temp_df4.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"
temp_df5.loc[temp_df5.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"
temp_df6.loc[temp_df6.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"
temp_df7.loc[temp_df7.value == 'SpliceSiteFinder, NNSPLICE', "value"] = "SpliceSiteFinder-like and NNSPLICE"


# %%
df_plot = temp_df.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?',index = "value", values = 0)
df_plot2 = temp_df2.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?', index = "value",values = 0)
df_plot3= temp_df3.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?',index = "value", values = 0)
df_plot4 = temp_df4.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?',index = "value", values = 0)
df_plot5 = temp_df5.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?',index = "value", values = 0)
df_plot6 = temp_df6.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?', index = "value",values = 0)
df_plot7 = temp_df7.groupby(['Are you completing this as a scientist, clinician, academic or other?',"value"]).size().reset_index().pivot(columns = 'Are you completing this as a scientist, clinician, academic or other?',index = "value", values = 0)

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05
    
    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    

plt.title(predictors_df.columns[1])
plt.ylabel('In-silico predictor')
plt.show()

# %%
import textwrap
plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot2.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    
                    
ax.set_yticklabels([textwrap.fill(e,16) for e in df_plot2.iloc[:,0].index.unique()])

plt.title(predictors_df.columns[2])
plt.ylabel('In-silico predictor')
plt.show()

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot3.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    

ax.set_yticklabels([textwrap.fill(e,16) for e in df_plot3.iloc[:,0].index.unique()])

plt.title(predictors_df.columns[3])
plt.ylabel('In-silico predictor')
plt.show()

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot4.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    

ax.set_yticklabels([textwrap.fill(e,16) for e in df_plot4.iloc[:,0].index.unique()])
plt.title(predictors_df.columns[4])
plt.ylabel('In-silico predictor')
plt.show()

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot5.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    
                
ax.set_yticklabels([textwrap.fill(e,16) for e in df_plot5.iloc[:,0].index.unique()])
plt.title(predictors_df.columns[5])
plt.ylabel('In-silico predictor')
plt.show()

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot6.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05
    
    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')    

ax.set_yticklabels([textwrap.fill(e,16) for e in df_plot6.iloc[:,0].index.unique()])
plt.title(predictors_df.columns[6])
plt.ylabel('In-silico predictor')
plt.show()

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot7.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')   

ax.set_yticklabels([textwrap.fill(e,16) for e in df_plot7.iloc[:,0].index.unique()])
plt.title(predictors_df.columns[7])
plt.ylabel('In-silico predictor')
plt.show()

# %%


# %%
import seaborn as sns
fig, ax = plt.subplots(3,2)
fig.set_figheight(15)
fig.set_figwidth(35)




sns.countplot(y = temp_df2.value, ax= ax[0][0], palette="Paired").set_title(predictors_df.columns[1])
sns.countplot(y = temp_df3.value, ax= ax[0][1], palette="Paired").set_title(predictors_df.columns[2])
sns.countplot(y = temp_df4.value, ax= ax[1][0], palette="Paired").set_title(predictors_df.columns[3])
sns.countplot(y = temp_df5.value, ax= ax[1][1], palette="Paired").set_title(predictors_df.columns[4])
sns.countplot(y = temp_df6.value, ax= ax[2][0], palette="Paired").set_title(predictors_df.columns[5])
sns.countplot(y = temp_df.value, ax= ax[2][1], palette="Paired").set_title(predictors_df.columns[0])

fig.tight_layout()
plt.show()

# %%
plt.figure(figsize = (7,5))
sns.countplot(y = temp_df5.value, palette="Paired").set_title(predictors_df.columns[5])
plt.ylabel('In-silico predictor')

# %%
temp_df_genes =df.iloc[:,27].str.split(";", 5, expand = True)
temp_df_genes["Job"] = df['Are you completing this as a scientist, clinician, academic or other?']

# %%
temp_df_genes =temp_df_genes.melt(id_vars = "Job")

# %%
temp_df_genes = temp_df_genes.drop(temp_df_genes.loc[(temp_df_genes.value == "None") | (temp_df_genes.value == "")].index)

# %%

df_plot = temp_df_genes.groupby(['Job', 'value']).size().reset_index().pivot(columns = 'Job', index = 'value', values = 0)

# %%

plt.rcParams["figure.figsize"] = (15,12)
total = df.shape[0]
ax = df_plot.plot(kind = 'barh', stacked = True)
for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy()
    x = x + width
    y = y + height / 2 - 0.05

    if width != 0:
        x, y = p.get_xy()
        if y > 0:
            ax.text(x + width/2,
                    y + height/2,
                    '{:.0f} %'.format(100*width/(total)),
                    horizontalalignment = 'center',
                    verticalalignment = 'center')   

plt.ylabel("Gene Level Metrics")
plt.show()

# %%

fig, ax = plt.subplots(figsize=(7, 5), dpi=96)
sns.countplot(y=df.iloc[:,5], hue = df.iloc[:,9], hue_order=["Yes", "No"],
            palette={"Yes": "Grey", "No": "Black"}, capsize=.1, ax=ax)
for bar_group, desaturate_value in zip(ax.containers, [0.5, 1]):
    for bar, color in zip(bar_group, ['blue', 'green', 'red', "orange"]):
        bar.set_facecolor(sns.desaturate(color, desaturate_value))
plt.show()

# %%
plt.figure(figsize = (10,7))
sns.countplot(y = temp_df_genes.value, palette="Paired").set_title("Gene Level Metrics")

# %%
plt.figure(figsize = (15,7))
sns.countplot(y = temp_df_genes.value, hue = temp_df_genes.Job, palette="Paired").set_title("Gene Level Metrics vs Profession")

# %%
plt.figure(figsize = (15,7))
sns.countplot(x = df['Would you be interested in using another automated system than the one you currently access?'], palette="Blues")
plt.xlabel('Interest in using an alternative automated classfication system')
plt.tight_layout()

# %%
plt.figure(figsize = (21,7))
sns.countplot(hue = df['Are you completing this as a scientist, clinician, academic or other?'], y = df['Would you be interested in using another automated system than the one you currently access?'], palette="Paired")
plt.ylabel('Interest in using an alternative automated classfication system')
plt.tight_layout()


# %%



