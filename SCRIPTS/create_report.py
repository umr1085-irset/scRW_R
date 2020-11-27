#############################################################
# Snakefile for single cell analysis workflow
# P. RIVAUD
# 2020/04
#############################################################

# In Python, Snakefile parameters can be accessed with the following:
# print(snakemake.input[1])
# print(snakemake.output["rds_sce"])

#########
# Imports
#########
from jinja2 import Environment, FileSystemLoader
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import json
from bokeh.palettes import all_palettes
import plotly.express as px
import math
import os

##############
# Get template
##############
env = Environment(loader=FileSystemLoader(snakemake.params['html_template_dir'])) # route to templates folder
template = env.get_template('template.html') # route to html template

#######################
# Render with variables
#######################
kwargs_dict = dict() # dictionary to store keyword arguments to render template
kwargs_dict['css']=snakemake.params['css_file']
kwargs_dict['plotly_js']=snakemake.input['plotlyjs_file']
kwargs_dict['report_icon']=os.path.join(snakemake.params['html_template_dir'],'icons/report_hex.svg')

##################################
# Data information section (step1)
##################################
kwargs_dict['project_name']=snakemake.params['project_name']
kwargs_dict['species']=snakemake.params['species']
kwargs_dict['number_of_samples'] = snakemake.params['nsamples']
kwargs_dict['step1_number_redundant_genes_removed'] = np.loadtxt(snakemake.input['datalog_step1_genes_rm'])

df = pd.read_csv(snakemake.input['datalog_df_step1'], header=0, index_col=0)

kwargs_dict['step1_total_cells'] = df.shape[0]
kwargs_dict['step1_avg_counts'] = round(df['sum'].mean(), 3)
kwargs_dict['step1_avg_detected'] = round(df['detected'].mean(), 3)

samples,ncells=[],[]
for idx,x in df.Sample.value_counts().iteritems():
	samples.append(idx)
	ncells.append(x)
fig = go.Figure(go.Bar(
            x=ncells,
            y=samples,
            orientation='h'))
fig.update_layout(
	title_text='Cells per sample',
	title_x=0.5,
	xaxis_title='Number of cells',
	#yaxis_title='Samples',
	height=110*snakemake.params['nsamples']
)
graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
kwargs_dict['step1_samples_cells_chart'] = graphJSON

# counts per cell
fig = go.Figure()
fig.add_trace(
    go.Histogram(x=df['sum'])
)
fig.add_shape(
        go.layout.Shape(type='line', xref='x', yref='paper', x0=kwargs_dict['step1_avg_counts'], y0=0, x1=kwargs_dict['step1_avg_counts'], y1=1, line={'dash': 'dash'}, name='Average counts'),
)
fig.update_layout(
	title_text='Counts per cell',
	title_x=0.5,
	xaxis_title='Counts',
	yaxis_title='Number of cells'
)
graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
kwargs_dict['step1_avg_counts_chart'] = graphJSON

# detected genes
fig = go.Figure()
fig.add_trace(
    go.Histogram(x=df['detected'])
)
fig.add_shape(
        go.layout.Shape(type='line', xref='x', yref='paper', x0=kwargs_dict['step1_avg_detected'], y0=0, x1=kwargs_dict['step1_avg_detected'], y1=1, line={'dash': 'dash'}, name='Average counts'),
)
fig.update_layout(
	title_text='Detected genes per cell',
	title_x=0.5,
	xaxis_title='Number of genes',
	yaxis_title='Number of cells'
)
graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
kwargs_dict['step1_genes_detected_chart'] = graphJSON

#######################
# Cell outliers (step2)
#######################
df = pd.read_csv(snakemake.input['datalog_df_step2'], header=0, index_col=0)
tmp = df['outlier'].value_counts()
kwargs_dict['step2_noutliers'] = tmp.loc[True]
kwargs_dict['step2_nremaining'] = tmp.loc[False]

fig = go.Figure()
colorsIdx = {True: 'crimson', False: '#636dfa'}
for outliervalue in [True, False]:
	dfSubset = df.loc[df['outlier'] == outliervalue]
	color = colorsIdx[outliervalue]
	fig.add_trace(
		go.Scattergl(
			x=dfSubset['PC1'],
			y=dfSubset['PC2'],
			mode='markers',
			name='%s (%d)' % (outliervalue,df['outlier'].value_counts()[outliervalue]),
			text=dfSubset.index,
			marker=dict(
				color=color,
				opacity=0.7
			)
		)
	)

fig.update_layout(
	title_text='Cells in PCA space computed from proportions of ribosomal and mitochondrial genes',
	title_x=0.5,
	xaxis_title='PC1',
	yaxis_title='PC2',
	legend_title='Outlier'
)
fig['layout']['yaxis']['scaleanchor']='x' # to unstretch the plot
graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
kwargs_dict['step2_cell_outliers_chart'] = graphJSON

########################
# Gene filtering (step4)
########################
df = pd.read_csv(snakemake.input['datalog_df_step4'], header=0, index_col=0)
df = df.sort_values('AveCount',ascending=False)
ntop = 25
genes = df.index.values[:ntop]
avg = df['AveCount'].values[:ntop]
fig = go.Figure(go.Bar(
            x=avg,
            y=genes,
            orientation='h'))
fig.update_layout(
	title_text='Top %d genes by average counts' % ntop,
	title_x=0.5,
	xaxis_title='Average counts',
	yaxis_title='Genes',
	height=600
)
graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
kwargs_dict['step4_avg_counts_chart'] = graphJSON

########################
# Doublet Finder (step5)
########################
samplelist = snakemake.params['samples']
subplot_titles= samplelist+['All samples']
ncols = 3
nrows = math.ceil((len(samplelist)+1)/ncols)
fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=subplot_titles)
for i, sample in enumerate(samplelist):
	with open('%s/%s_barcodes_singlets.txt' % (snakemake.params['doubletfinderdir'],sample)) as f:
		for nsinglets, l in enumerate(f,1):
			pass
	with open('%s/%s_barcodes_doublets.txt' % (snakemake.params['doubletfinderdir'], sample)) as f:
		for ndoublets, l in enumerate(f,1):
			pass
	labels=['Singlets', 'Doublets']
	fig.add_trace(
		go.Bar(x=labels, y=[nsinglets, ndoublets], marker_color=['#636dfa','crimson']),
		row=int(i/3)+1, col=(i%3)+1
	)

i+=1
with open('%s/ALL_barcodes_singlets.txt' % snakemake.params['doubletfinderdir']) as f:
		for nsinglets, l in enumerate(f,1):
			pass
ndoublets = kwargs_dict['step2_nremaining']-nsinglets
fig.add_trace(
	go.Bar(x=labels, y=[nsinglets, ndoublets], marker_color=['#636dfa','crimson']),
	row=int(i/3)+1, col=(i%3)+1
)
fig.update_layout(
	title_text="Singlets and doublets amounts", 
	title_x=0.5,
	showlegend=False,
	height=450*nrows
)
graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
kwargs_dict['step5_doubletfinder_chart'] = graphJSON
		
###########################################
# Normalizations (step7) SEURAT_SCTRANSFORM
###########################################
def create_umap_fig(df, col, norm):
	fig = go.Figure()
	for value in df[col].unique():
		dfSubset = df.loc[df[col] == value]
		fig.add_trace(
			go.Scattergl(
				x=dfSubset['UMAP_1'],
				y=dfSubset['UMAP_2'],
				mode='markers',
				name='%s' % value,
				text=dfSubset.index,
				marker=dict(
					opacity=0.7
				)
			)
		)
	
	if col == 'Sample':
		title = "UMAP colored by sample, from %s normalization" % norm.replace('_',' ').title()
	elif col == 'Phase':
		title = "UMAP colored by cell cycle phase, from %s normalization" % norm.replace('_',' ').title()
	elif col == 'clust':
		title = "UMAP colored by cluster, from %s normalization" % norm.replace('_',' ').title()

	fig.update_layout(
		title_text=title,
		title_x=0.5,
		showlegend=True,
		height=800,
		width=900
	)

	return fig

for norm in ['seurat_sctransform', 'seurat_lognorm', 'scran_deconvolution']:
	try:
		datalogfile = '%sdatalogstep7_%s_df.csv' % (snakemake.params['basedir_step7'],norm)
		df = pd.read_csv(datalogfile, header=0, index_col=0)
		s = '<h1>Normalization: %s</h1>' % norm.replace('_',' ').title()
		kwargs_dict['div_%s' % norm] = s

		fig = create_umap_fig(df, 'Sample', norm)
		graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
		kwargs_dict['step7_%s_chart_sample' % norm] = graphJSON
		
		fig = create_umap_fig(df, 'Phase', norm)
		graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
		kwargs_dict['step7_%s_chart_phase' % norm] = graphJSON
		
		fig = create_umap_fig(df, 'clust', norm)
		graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
		kwargs_dict['step7_%s_chart_cluster' % norm] = graphJSON
			
	except:
		pass

# Render HTML template
output_from_parsed_template = template.render(**kwargs_dict)

# to save the results
with open(snakemake.output['html_report'], "w") as fout:
    fout.write(output_from_parsed_template)

with open(snakemake.output['step_complete'], 'w') as fp: 
    pass

'''
<table class="table table-hover">
'''