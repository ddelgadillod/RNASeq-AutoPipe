## Librerias Shiny app
import pandas as pd
from shiny.types import ImgData
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
## Librerias flujo

import os
import pickle as pkl
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import pandas as pd
import numpy as np
import scanpy as sc
import requests
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage


## defiiendo funciones

def get_gene_symbol(ensembl_gene_id):
    # URL de la API REST de Ensembl
    url = f'https://rest.ensembl.org/lookup/id/{ensembl_gene_id}?content-type=application/json'
    response = requests.get(url).json()
    print(response)
    # Si la respuesta contiene el símbolo, lo devuelve, sino devuelve 'N/A'
    return response.get('display_name', 'N/A')

def categorize(row):
    if row["padj"] < 0.05 and row["log2FoldChange"] > 1:
        return "Upregulated"
    elif row["padj"] < 0.05 and row["log2FoldChange"] < -1:
        return "Downregulated"
    else:
        return "Non-significant"


## Parametros definidos por el usuario

SAVE = True  # whether to save the outputs of this notebook
countsPath = "/home/alejo/Avena-PyDS2/counts-8000.csv"
metadataPath = "/home/alejo/Avena-PyDS2/metadata.csv"


sampleIDCol = 'shortName'
metaCols = ['shortName', 'meta1', 'meta2', 'meta3','replicate']

condition = 'meta1'

minCounts = 20

### Contrast, construir a partir de condition y levels de condition
contrastList = ['meta1', 'Estres', 'Control']

### Coefficient string construir a partir de meta1 y condicion contrastante 
coeffStr = "meta1[T.Estres]"

### Samples a excluir

Samples2Exc = []

minBaseMean = 10

Geneid = "name"

N_top_genes = 20


min_expresed_genes = 10


## Importando datos

# counts = pd.read_csv(countsPath)
# counts = counts.set_index(Geneid)
# counts = counts[counts.sum(axis = 1) > 0]
# counts = counts.T

#counts

# metadata = pd.read_csv(metadataPath)
# metadata = metadata[metaCols]
# metadata = metadata.set_index(sampleIDCol)


# samples_to_keep = ~metadata[condition].isna()
# counts_df = counts.loc[samples_to_keep]

# metadata = metadata.loc[samples_to_keep]

# genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= minCounts]

# counts_df = counts_df[genes_to_keep]


## DESeq2

inference = DefaultInference(n_cpus=8)

global_data = {
    'counts': None,
    'metadata': None,
    'metadata_raw': None,
    'counts_df': None,
    'dds': None,
    'dds': None,
    'LFC' : None,
    'ds' : None,
    'res': None,
    'sigs': None,
    'normalized_counts' : None,
    'top_gene_counts' : None,
    'log1p_counts' : None,
    'top_genes' : None,
    'vst' : None,
    'vst_distance' : None,
    'log_counts' : None,
    'log_counts_long' : None,
    'correlation_matrix' : None,
    'expressed_genes' : None,
    'padj_thresh' : None,
    'lfc_null' : None,
    'colorsMA' : None,
    'disps' : None,
    'colorsDisp' : None,
    'legend_labels' : None,
    'pca' : None,
    'pca_result' : None,
    'pca_df' : None,
    'dds_sigs' : None,
    'grapher' : None,
    'resultsF' : None,
    'top_upregulated' : None,
    'top_downregulated' : None
}



# Define la interfaz de usuario
app_ui = ui.page_fluid(
    ui.h2("Flujo de análisis de datos RNA Seq"),
    ui.input_file("count_upload", "Cargar archivo CSV", accept=[".csv"]),
    ui.input_action_button("btn_showcounts", "Visualizar Counts"),
    ui.output_text("message_output"),
    ui.output_data_frame("raw_counts"),
    ui.input_file("meta_upload", "Cargar archivo Metadatos", accept=[".csv"]),
 
    ui.input_action_button("btn_showmeta", "Visualizar Metadatos"),
    ui.output_text("message_output_meta"),
    ui.output_data_frame("raw_metadata"),
    ui.input_action_button("btn_calc", "Calcular", disabled = False),
    
    ui.output_data_frame("results_df"), #height="400px",
    
    ui.input_action_button("btn_plot", "Graficar", disabled = True),
    ## Heatmap best genes * muestra
    ui.output_plot("HM_bestSample", width='60%', height='600px' ),
    ui.output_plot("HM_distancias", width='60%', height='800px' ),

    ui.output_plot("BP_conteos", width='60%', height='600px' ),
    ui.output_plot("VP_conteos", width='60%', height='600px' ),

    ui.output_plot("HM_correlaciones", width='60%', height='800px' ),

    ui.output_plot("BP_expresion", width='60%', height='600px' ),

    ui.output_plot("MA", width='60%', height='600px' ),

    ui.output_plot("ModDispersion", width='60%', height='600px' ),

    ui.output_plot("PCA1", width='60%', height='600px' ),
    ui.output_plot("PCA2", width='60%', height='600px' ),
    ui.output_plot("PCA3", width='60%', height='600px' ),

    ui.output_plot("HMmain", width='70%', height='800px' ),

    ui.output_plot("Vulcano", width='70%', height='700px' ),


    #ui.output_image("MA_png"),
    #ui.output_table("results_df"),
)

# Define el servidor
def server(input, output, session):
    data1 = reactive.Value(None)
    data2 = reactive.Value(None)


    @reactive.Effect
    @reactive.event(input.btn_showcounts)
    def load_file():
        # Obtener archivo subido
        uploaded_file = input.count_upload()
        if not uploaded_file:
            return  # No hacer nada si no hay archivo cargado

        # Leer el archivo CSV
        try:
            counts = pd.read_csv(uploaded_file[0]["datapath"])
            counts = counts.set_index(Geneid)
            counts = counts[counts.sum(axis = 1) > 0]
            counts = counts.T
            data1.set(counts.T.head(10))
            global_data["counts"] = counts
        except Exception as e:
            data1.set(None)
            print(f"Error al cargar el archivo: {e}")

    @output
    @render.text
    def message_output():
        if data1() is None:
            return "Por favor, carga un archivo CSV y presiona el botón para mostrar el contenido."
        return ""


    #@output
    @render.data_frame
    def raw_counts():
        if data1() is not None:
            return data1()
        return None


    @reactive.Effect
    @reactive.event(input.btn_showmeta)
    def load_file():
        # Obtener archivo subido
        uploaded_file = input.meta_upload()
        if not uploaded_file:
            return  # No hacer nada si no hay archivo cargado

        # Leer el archivo CSV
        try:
            metadata = pd.read_csv(uploaded_file[0]["datapath"])
            global_data["metadata_raw"] = metadata
            data2.set(metadata.head(10))
            metadata = metadata[metaCols]
            metadata = metadata.set_index(sampleIDCol)
            counts = global_data["counts"]
            samples_to_keep = ~metadata[condition].isna()
            counts_df = counts.loc[samples_to_keep]

            metadata = metadata.loc[samples_to_keep]

            genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= minCounts]

            counts_df = counts_df[genes_to_keep]
            
            
            global_data["metadata"] = metadata
            global_data["counts_df"] = counts_df
        except Exception as e:
            data2.set(None)
            print(f"Error al cargar el archivo: {e}")

    @output
    @render.text
    def message_output_meta():
        if data2() is None:
            return "Por favor, carga un archivo CSV y presiona el botón para mostrar el contenido."
        return ""


    #@output
    @render.data_frame
    def raw_metadata():
        if data2() is not None:
            return data2()
        return None   


    @reactive.effect
    @reactive.event(input.btn_calc)
    def enable_graficar():
        if input.btn_calc:
            ui.update_action_button("btn_plot", disabled=False)
        else:
            ui.update_action_button("btn_plot", disabled=True)



    @render.data_frame
    @reactive.event(input.btn_calc)
    def results_df():
        counts_df = global_data['counts_df']
        metadata = global_data['metadata']
        inference = DefaultInference(n_cpus=8)
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design="~"+ condition,
            refit_cooks=True,
            inference=inference,
            # n_cpus=8, # n_cpus can be specified here or in the inference object
        )

        dds.deseq2()

        LFC = dds.varm["LFC"]



        ds = DeseqStats(dds, contrast=contrastList, inference=inference)

        ds.summary(lfc_null=0.1, alt_hypothesis="greaterAbs")

        ds.lfc_shrink(coeff=coeffStr)

        print("Results to dataframe")
        res = ds.results_df
        res = res[res.baseMean >= minBaseMean]
        sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]



        normalized_counts = pd.DataFrame(dds.layers['normed_counts'])
        normalized_counts.columns = counts_df.columns
        normalized_counts.index = counts_df.index
        print(normalized_counts.head())

        print(dds)


        print(res.head())

        dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
        dds.layers['log1p']


        log1p_counts = pd.DataFrame(dds.layers['log1p'])
        log1p_counts.columns = counts_df.columns
        log1p_counts.index = counts_df.index
        print(log1p_counts)

        top_genes = log1p_counts.T.mean(axis=1).sort_values(ascending=False).head(N_top_genes).index
        top_gene_counts = log1p_counts.T.loc[top_genes]
        print(top_genes)


        print(top_gene_counts)

        top_gene_counts = top_gene_counts[metadata.index]
        #global_data['top_gene_counts'] = top_gene_counts
        print(top_gene_counts)

        dds.vst()

        vst = pd.DataFrame(dds.layers['vst_counts'])
        vst.columns = counts_df.columns
        vst.index = counts_df.index
        vst = vst.T
        print(vst)

        vst_distance = pd.DataFrame(
            squareform(pdist(vst.T, metric="euclidean")),  # Calcular distancias entre muestras
            index=vst.columns,  # Etiquetas de las muestras
            columns=vst.columns
        )
        vst_distance

        print(vst_distance)

        log_counts = np.log2(normalized_counts + 1)
        print(log_counts)

        log_counts_long = log_counts.T.melt(var_name="Muestra", value_name="Log2(Conteos + 1)")
        print(log_counts_long)

        correlation_matrix = normalized_counts.T.corr()
        print(correlation_matrix)

        expressed_genes = (normalized_counts.T > min_expresed_genes).sum(axis=0)
        print(expressed_genes)

        padj_thresh = 0.05
        lfc_null = 0
        colorsMA = res["padj"].apply(lambda x: "darkred" if x < padj_thresh else "gray")
        print(colorsMA)

        disps = [
            dds.varm["genewise_dispersions"],
            dds.varm["dispersions"],
            dds.varm["fitted_dispersions"],
        ]
        print(disps)

        if len(disps) == 3:
            colorsDisp = "kbr"
        else:
            colorsDisp = "kbrcmyg"

        legend_labels = ["Estimated", "Final", "Fitted"]

        pca = PCA(n_components=3)
        pca_result = pca.fit_transform(normalized_counts)
        print(pca_result)

        pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2", "PC3"])
        pca_df["condition"] = metadata[condition].values  # Agregar la condición experimental
        pca_df["Sample"] = metadata.index
        print(pca)


        dds_sigs = dds[:, sigs.index]
        print(dds_sigs)

        grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                            index=dds_sigs.var_names, columns=dds_sigs.obs_names)

        print(grapher)

        resultsF = res.dropna(subset=["log2FoldChange", "padj"])
        resultsF["significance"] = (resultsF["padj"] < 0.05) & (abs(resultsF["log2FoldChange"]) > 1)
        print(resultsF)

        resultsF["category"] = resultsF.apply(categorize, axis=1)

        top_upregulated = resultsF[resultsF["category"] == "Upregulated"].sort_values(by="log2FoldChange", ascending=False).head(5)
        top_downregulated = resultsF[resultsF["category"] == "Downregulated"].sort_values(by="log2FoldChange").head(5)

        global_data['dds'] = dds
        global_data['LFC'] = LFC
        global_data['ds'] = ds
        global_data['res'] = res
        global_data['sigs'] = sigs
        global_data['normalized_counts'] = normalized_counts
        global_data['top_gene_counts'] = top_gene_counts
        global_data['log1p_counts'] = log1p_counts
        global_data['top_genes'] = top_genes
        global_data['vst'] = vst
        global_data['vst_distance'] = vst_distance
        global_data['log_counts'] = log_counts
        global_data['log_counts_long'] = log_counts_long
        global_data['correlation_matrix'] = correlation_matrix
        global_data['expressed_genes'] = expressed_genes
        global_data['padj_thresh'] = padj_thresh
        global_data['lfc_null'] = lfc_null
        global_data['colorsMA'] = colorsMA
        global_data['disps'] = disps
        global_data['colorsDisp'] = colorsDisp
        global_data['legend_labels'] = legend_labels
        global_data['pca'] = pca
        global_data['pca_result'] = pca_result
        global_data['pca_df'] = pca_df
        global_data['dds_sigs'] = dds_sigs
        global_data['grapher'] = grapher
        global_data['resultsF'] = resultsF
        global_data['top_upregulated'] = top_upregulated
        global_data['top_downregulated'] = top_downregulated


        return render.DataGrid(res.head(10))
    

    @render.plot(alt="Heatmap") 
    @reactive.event(input.btn_plot)
    def HM_bestSample():
        metadata = global_data['metadata']
        top_gene_counts = global_data['top_gene_counts']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 8))

        ax = sns.heatmap(
            top_gene_counts, 
            cmap="RdYlBu_r", 
            cbar=True,
            xticklabels=metadata.index,  # Usar "shortName" (índice de design_df) como etiquetas de columnas
            yticklabels=False,           # Ocultar nombres de filas (genes)
            linewidths=0.5
            )

        plt.title("Heatmap de los " + str(N_top_genes) + " genes más expresados")
        plt.xlabel("Muestras")
        plt.ylabel("Genes")

        return ax 
    
    @render.plot(alt="Heatmap")
    @reactive.event(input.btn_plot)
    def HM_distancias():
        vst_distance = global_data['vst_distance']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 8))   
        ax = sns.clustermap(
            vst_distance,
            metric="euclidean",       # Distancia usada para agrupar
            method="ward",            # Método de agrupamiento jerárquico
            cmap="RdYlBu_r",           # Esquema de colores
            figsize=(10, 10),         # Tamaño del gráfico
            xticklabels=True,         # Mostrar etiquetas de columnas
            yticklabels=True,
            fmt=".2f", # Mostrar etiquetas de filas
        )

        # Título del gráfico
        #plt.title("Heatmap de Distancias varianza estabilizada")
        #plt.show()
        return ax

    @render.plot(alt="BoxPlot")
    @reactive.event(input.btn_plot)
    def BP_conteos():
        log_counts = global_data['log_counts']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(data=log_counts.T, palette="Set3", showfliers=False)
        plt.title("Distribución de los conteos normalizados por muestra")
        plt.xlabel("Muestras")
        plt.ylabel("Log2(Conteos Normalizados + 1)")
        plt.xticks(rotation=90)
        #plt.show()

        return ax

    @render.plot(alt="ViolinPlot")
    @reactive.event(input.btn_plot)
    def VP_conteos():
        log_counts_long = global_data['log_counts_long']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.violinplot(
            x="Muestra",
            y="Log2(Conteos + 1)",
            data=log_counts_long,
            palette="Set3"
        )
        plt.title("Distribución de los conteos normalizados (Violin Plot)")
        plt.xlabel("Muestras")
        plt.ylabel("Log2(Conteos Normalizados + 1)")
        plt.xticks(rotation=90)
                #plt.show()

        return ax

    @render.plot(alt="Heatmap")
    @reactive.event(input.btn_plot)  
    def HM_correlaciones():
        correlation_matrix = global_data['correlation_matrix']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(
            correlation_matrix, 
            annot=True, 
            fmt=".2f", 
            cmap="RdYlBu_r", 
            square=True
        )
        plt.title("Matriz de correlaciones entre muestras")

        return ax 

    @render.plot(alt="BarPlot")
    @reactive.event(input.btn_plot)  
    def BP_expresion():
        expressed_genes = global_data['expressed_genes']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.barplot(x=expressed_genes.index, y=expressed_genes.values, palette="Set3")
        plt.title("Número de genes expresados por muestra (> " + str(min_expresed_genes) + "reads)")
        plt.xlabel("Muestras")
        plt.ylabel("Número de genes expresados")
        plt.xticks(rotation=90)


        return ax 

    @render.plot(alt="MA-Plot")
    @reactive.event(input.btn_plot)  
    def MA():
        res = global_data['res']
        colorsMA = global_data['colorsMA']
        padj_thresh = global_data['padj_thresh']
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = plt.scatter(
            x=res["baseMean"],
            y=res["log2FoldChange"],
            c=colorsMA,
            alpha = 0.5,
            s = 20
        )
        plt.xscale("log")
        plt.xlabel("Promedio de conteos normalizados")
        plt.ylabel("log2 fold change")
        plt.axhline(padj_thresh, color="red", alpha=0.5, linestyle="--", zorder=3)
        plt.axhline(-padj_thresh, color="red", alpha=0.5, linestyle="--", zorder=3)
        plt.tight_layout()
        return ax 

    @render.plot(alt="MA-Plot")
    @reactive.event(input.btn_plot)  
    def ModDispersion():
        dds = global_data['dds']
        disps = global_data['disps']
        colorsDisp = global_data['colorsDisp']
        legend_labels = global_data['legend_labels']

        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10, 6))


        for disp, color in list(zip(disps, colorsDisp)):
            plt.scatter(x=dds.varm["_normed_means"], y=disp, c=color, alpha = 0.5, s = 0.6)
        plt.yscale("log")
        plt.xscale("log")
        # label legend + axes
        plt.legend(legend_labels, loc="best")
        plt.xlabel("Promedio de conteos normalizados")
        plt.ylabel("Dispersión")
        plt.tight_layout()

        return 

    @render.plot(alt="PCA1")
    @reactive.event(input.btn_plot)  
    def PCA1():
        pca = global_data['pca']
        pca_df = global_data['pca_df']

        sns.set(style="whitegrid")
        plt.figure(figsize=(8, 6))

        ax = sns.scatterplot(
            data=pca_df,
            x="PC1",
            y="PC2",
            hue="condition",  # Cambia 'condition' según tu variable de interés
            palette="Set3",
            s=100
        )

        for i, label in enumerate(pca_df['Sample']):
            plt.text(
                pca_df['PC1'][i] + 0.02,  # Desplazamiento en X
                pca_df['PC2'][i] + 0.02,  # Desplazamiento en Y
                label,
                fontsize=9
            )

        plt.title("PCA de datos normalizados")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} varianza explicada)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} varianza explicada)")
        plt.legend(title="Condición")
        return ax

 

    @render.plot(alt="PCA2")
    @reactive.event(input.btn_plot)  
    def PCA2():
        pca = global_data['pca']
        pca_df = global_data['pca_df']
        sns.set(style="whitegrid")
        plt.figure(figsize=(8, 6))

        ax = sns.scatterplot(
            data=pca_df,
            x="PC1",
            y="PC3",
            hue="condition",  # Cambia 'condition' según tu variable de interés
            palette="Set3",
            s=100
        )
        
        for i, label in enumerate(pca_df['Sample']):
            plt.text(
                pca_df['PC1'][i] + 0.02,  # Desplazamiento en X
                pca_df['PC3'][i] + 0.02,  # Desplazamiento en Y
                label,
                fontsize=9
            )

        plt.title("PCA de datos normalizados")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} varianza explicada)")
        plt.ylabel(f"PC3 ({pca.explained_variance_ratio_[2]:.2%} varianza explicada)")
        plt.legend(title="Condición")
        return ax

    @render.plot(alt="PCA3")
    @reactive.event(input.btn_plot)  
    def PCA3():
        pca = global_data['pca']
        pca_df = global_data['pca_df']
        sns.set(style="whitegrid")
        plt.figure(figsize=(8, 6))

        ax = sns.scatterplot(
            data=pca_df,
            x="PC2",
            y="PC3",
            hue="condition",  # Cambia 'condition' según tu variable de interés
            palette="Set3",
            s=100
        )
        
        for i, label in enumerate(pca_df['Sample']):
            plt.text(
                pca_df['PC2'][i] + 0.02,  # Desplazamiento en X
                pca_df['PC3'][i] + 0.02,  # Desplazamiento en Y
                label,
                fontsize=9
            )

        plt.title("PCA de datos normalizados")
        plt.xlabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} varianza explicada)")
        plt.ylabel(f"PC3 ({pca.explained_variance_ratio_[2]:.2%} varianza explicada)")
        plt.legend(title="Condición")
        return ax

    @render.plot(alt="main Heatmap")
    @reactive.event(input.btn_plot)  
    def HMmain():
        grapher = global_data['grapher']
        sns.set(style="whitegrid")
        plt.figure(figsize=(20, 15))
        ax = sns.clustermap(grapher, z_score=0, cmap = 'RdYlBu_r')
        return ax

    @render.plot(alt="Vulcano")
    @reactive.event(input.btn_plot)  
    def Vulcano():
        grapher = global_data['grapher']
        resultsF = global_data['resultsF']
        top_upregulated = global_data['top_upregulated']
        top_downregulated = global_data['top_downregulated']

        plt.figure(figsize=(12, 8))
        ax = sns.scatterplot(
            data=resultsF,
            x="log2FoldChange",
            y=-np.log10(resultsF["padj"]),
            hue="category",  # Colorear según la categoría
            palette={"Upregulated": "red", "Downregulated": "blue", "Non-significant": "grey"},
            legend="full",
            s=80
        )

        # Añadir líneas de referencia
        plt.axhline(-np.log10(0.05), color="black", linestyle="--", label="p-adj=0.05")
        plt.axvline(-1, color="black", linestyle="--", label="log2FoldChange=-1")
        plt.axvline(1, color="black", linestyle="--", label="log2FoldChange=1")

        # Etiquetar los 5 más upregulated
        for i, row in top_upregulated.iterrows():
            plt.text(
                x=row["log2FoldChange"], 
                y=-np.log10(row["padj"]) + 0.2,  # Desplazamiento en Y
                s=i,  # Nombre del gen
                fontsize=9,
                color="darkred"
            )

        # Etiquetar los 5 más downregulated
        for i, row in top_downregulated.iterrows():
            plt.text(
                x=row["log2FoldChange"], 
                y=-np.log10(row["padj"]) + 0.2,  # Desplazamiento en Y
                s=i,  # Nombre del gen
                fontsize=9,
                color="darkblue"
            )

        # Etiquetas y título
        plt.title("Volcano Plot: Upregulated, Downregulated, and No significativos", fontsize=16)
        plt.xlabel("Log2 Fold Change", fontsize=14)
        plt.ylabel("-Log10 p-value ajustado", fontsize=14)
        plt.legend(title="Category")
        return ax

# Crear la aplicación
app = App(app_ui, server)