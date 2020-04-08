"""
Creates a Dash application that provides 3D visualization of a KINC network.

This script accepts the following arguments:

    --net : (required) The path to the KINC-derived network file

    --gem : (retuired) The path to the log2 transformed Gene Expression Matrix.

    --amx : (required) The path to the tab-delimited annotation matrix.
            The matrix must have at least one column that contains a unique
            list of sample names.  The column must be named 'Sample'.

    --color_col : (required) The name of the column in the annotation matrix
                  The contains the categories that should be colored in the
                  display.

"""

import argparse
import numpy as np
import pandas as pd
import igraph as ig
import plotly as py
import seaborn as sns
import plotly.graph_objects as go
import networkx as nx
import dash
import dash_core_components as dcc
import dash_html_components as html
import os
import json
import re





def load_network(file_path):
    """
    Imports the KINC-generated network file (either full or Tidy versions).

    file_path : The path to the network file.

    return : A pandas dataframe containing the network.
    """

    net = pd.read_csv(file_path, sep="\t")
    return net





def load_gem(file_path):
    """
    Imports the tab-delimited Gene Expression Matrix (GEM).

    GEM files can be generated from RNA-seq data using GEMmaker.

    file_path : The path to the GEM file.  The file should be log2 transformed.

    return : A pandas dataframe containing the GEM.
    """

    gem = pd.read_csv(file_path, sep="\t")
    return gem





def load_amx(file_path, sample_col = 'Sample'):
    """
    Imports the tab-delimited annotation matrix (amx).

    The matrix must have at least one column that contains a unique list of
    sample names.

    file_path :  The path to the annotation matrix.

    sample_col : The name of the column that contains the sample names. Defaults
                 to 'Sample'

    return : A pandas dataframe containing the annotation matrix.
    """
    amx = pd.read_csv(file_path, sep="\t")
    amx.index = amx[sample_col]
    return amx





def get_iGraph(net):
    """
    Converts the KINC network dataframe into an iGraph object.

    Igraph objects are handly for performing network statistics such as
    transitivity and degree calculations.

    net :  The network dataframe created by the load_network function.

    return : An igraph object of the network loaded with the source, target and
             Similarity_Score (as the weight)
    """
    g = ig.Graph()

    # Add the nodes
    v = pd.concat([net['Source'], net['Target']]).unique()
    g.add_vertices(v)

    # Add the edges
    g.add_edges(net[['Source', 'Target']].values)

    # Add the edge w
    g.es['weight'] = net['Similarity_Score']

    return g





def calculate_2d_layout(net):
    """
    Calculates a typical 2D layout for the network.

    The first time this function is called on a network it may take some time
    depending on the size of the network.  The layout is saved in a file named
    'glayout.txt' in the current working directory. On subsequent runs of
    this program that file is imported if it exists.

    net :  The network dataframe created by the load_network function.

    return : a Pandas dataframe containing the layout coordinates for
             the nodes in the network. The dataframe contains X, and Y
             dimenstional coordinates.
    """

    g = get_iGraph(net)
    t = pd.DataFrame(g.transitivity_local_undirected())[0].values
    d = pd.DataFrame(g.degree(), index=g.vs['name'], columns=['Degree'])

    G = nx.Graph()
    net['Weight'] = np.abs(net['Similarity_Score'])
    G.add_weighted_edges_from(net[['Source','Target','Weight']].values)
    if (not os.path.exists('glayout.txt')):
        glayout = pd.DataFrame(nx.drawing.layout.kamada_kawai_layout(G)).transpose()
        glayout.columns = ['X', 'Y']
        glayout = pd.concat([glayout, d, t], axis=1, sort=False)
        glayout.to_csv('glayout.txt')
    else:
        glayout = pd.read_csv('glayout.txt', index_col=0)

    return glayout





def get_zlayer_bins(net):
    """
    Calculates a set of bins using the Similarity score.

    It is from these bins that the edges and nodes of the network will be
    stacked in the z-axis of the 3D plot.

    net :  The network dataframe created by the load_network function.

    return : a numpy array containing the list of bin values.
    """
    return np.around(np.abs(net['Similarity_Score']), decimals=2)





def get_vertex_zlayers(net, glayout):
    """
    Uses the 2D layout and calculates the Z-coordinate for the nodes.

    net :  The network dataframe created by the load_network function.

    glayout : The dataframe containing the 2D layout of the nodes.

    return : A Pandas dataframe containing the X, Y and Z coordinates for the
             nodes as well as the Degree and CC (clustering coefficient) for
             each node.
    """
    net['Bin'] = get_zlayer_bins(net)
    net['Weight'] = np.abs(net['Similarity_Score'])

    def find_vlayers(row, vtype='Source'):
        node = glayout.loc[row[vtype]]
        sbin = row['Bin']
        return(row[vtype], node['X'], node['Y'], sbin, node['Degree'], node['CC'])

    lsource = net.apply(find_vlayers, vtype='Source', axis=1)
    ltarget = net.apply(find_vlayers, vtype='Target', axis=1)

    vlayers = pd.DataFrame.from_records(lsource.append(ltarget).values, columns=['Vertex', 'X', 'Y', 'Z', 'Degree', 'CC'])
    vlayers.dropna(inplace=True)
    vlayers = vlayers[vlayers.duplicated() == False]
    vlayers.head()

    vlayers = vlayers.groupby(by=['Vertex']).apply(lambda g: g[g['Z'] == g['Z'].max()])
    return vlayers





def get_edge_zlayers(net, glayout):
    """
    Uses the 2D layout and calculates the Z-coordinate for the edges.

    Edges are drawn as lines in the 3D scatterplot, therefore this function
    calculates the start and stop coordinates for the edges in the format
    required by the scatter3d viewer.

    net :  The network dataframe created by the load_network function.

    glayout : The dataframe containing the 2D layout of the nodes.

    return : A Pandas dataframe containing the X, Y and Z coordinates arrays
             for the edges as well as Source, Target and Samples values from
             the original network.  The X, Y and Z coordiantes are tuples.
    """

    def place_elayers(row):
        sbin = row['Bin']
        source = glayout.loc[row["Source"]]
        target = glayout.loc[row["Target"]]
        return([[source['X'], target['X'], None],
                [source['Y'], target['Y'], None],
                [sbin, sbin, None],row["Source"],row["Target"],row["Samples"]])

    ledge = net.apply(place_elayers, axis=1)

    elayers = pd.DataFrame.from_records(ledge, columns=['X', 'Y', 'Z', 'Source', 'Target', 'Samples'])
    elayers.dropna(inplace=True)
    elayers['name'] = elayers['Source'] + " (co) " + elayers['Target']

    return elayers





def create_network_plot(net, vlayers, elayers):
    """
    Uses Plotly to create the interactive 3D visualization of the network.

    This function uses the Scatter3D plot to draw the network.  The axes are
    hidden so it appears as a typical network view.  It defaults to
    a straight on view as the network would be seen in a typical 2D viewer like
    Cytoscape.

    net :  The network dataframe created by the load_network function.

    vlayers : The dataframe containing the 3D coordinates for the nodes.

    elayers : The dataframe containing the 3D coordinates for the edges.

    return : a Plotly figure object.
    """

    fig1 = go.Figure(data=[go.Scatter3d(
                   x=vlayers['X'],
                   y=vlayers['Y'],
                   z=vlayers['Z'],
                   mode='markers',
                   marker=dict(symbol='circle', size=np.log10(vlayers['Degree'])*4,
                               color=vlayers['Z'], colorscale='Viridis'),
                   text="Node: " + vlayers['Vertex'],
                   hoverinfo='text')])

    # Calculate a color per layer
    net['Bin'] = get_zlayer_bins(net)
    bins = np.flip(np.sort(net['Bin'].unique()))
    layer_colors = pd.DataFrame({
        'bin': bins,
        'color' : np.flip(sns.color_palette('viridis', bins.size).as_hex())
    })

    # Add edge traces to the figure, one each per bin.
    for bin_val in bins:
        layer_indexes = (net['Bin'] == bin_val)
        layer_color = layer_colors[layer_colors['bin'] == bin_val]['color'].values[0]

        # Reformat the elayers for use by the Scatter3d function.
        eX = np.hstack(elayers['X'][layer_indexes])
        eY = np.hstack(elayers['Y'][layer_indexes])
        eZ = np.hstack(elayers['Z'][layer_indexes])

        # Calculate the array to represent colors.
        eZdf = pd.DataFrame.from_records(
                  elayers['Z'],
                  columns=['Source', 'Target', 'NULL'])[layer_indexes]
        eZdf['NULL'] = eZdf['Source']
        eC = np.hstack(eZdf.values)

        # Create the scatterplot containing the lines for edges.
        fig1.add_trace(go.Scatter3d(
                       x=eX,
                       y=eY,
                       z=eZ,
                       mode='lines',
                       line=dict(color=layer_color, width=1),
                       text="Edge: " + elayers["name"],
                       hoverinfo='text'))

    # ---------------------------------------------------------------
    #  Add a slider for the network viewer
    # ---------------------------------------------------------------

    steps = []
    num_layers = len(fig1.data)
    for i in range(num_layers):
        step = dict(
            method="restyle",
            # Disable all layers.
            args=["visible", [False] * num_layers if (i > 0) else [True] * num_layers],
            label='all' if (i == 0) else bins[i-1]
        )
        # Turn on the layers for this step. It should always be the first
        # layer and i'th layer.
        step["args"][1][0] = True
        for j in range(i):
            step["args"][1][j] = True

        # Set the label.
        steps.append(step)

    fig1.update_layout(
        height=600,
        title=dict(text = "3D Network View", font = dict(color='#FFFFFF')),
        showlegend=False,
        margin=dict(l=10, r=10, t=30, b=10),
        paper_bgcolor = "#000000",
        scene=dict(
          aspectmode="cube",
          xaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False),
          yaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=False, title='', showspikes=False),
          zaxis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False,
                     showticklabels=True, title='', showspikes=False, color="#FFFFFF")
        ),
        hovermode='closest',
        annotations=[dict(showarrow=False,
                        text="",
                        xref='paper', yref='paper',
                        x=0, y=0.1, xanchor='left', yanchor='bottom', font=dict(size=14))
                    ],
        sliders=[dict(
            active=0,
            currentvalue={"prefix": "Similarity: "},
            pad={"t": 50},
            steps=steps,
            font=dict(color = '#FFFFFF'),
            tickcolor='#FFFFFF',
            len=1)],
    )

    # We want an orthographic layout so that when looking above the edges line up
    # with the nodes.
    fig1.layout.scene.camera.projection.type = "orthographic"
    fig1.layout.scene.camera.eye = dict(x=0., y=0., z=2)

    return fig1





def create_expression_scatterplot(gem, amx, elayers, color_col, edge_index=0):
    """
    Uses Plotly to create the interactive 3D scatterplot of co-expression

    This function uses the Scatter3D plot to draw the co-expression scatterplot.
    It defaults to a straight on view but can be interactively rotated,
    panned, etc.

    net :  The network dataframe created by the load_network function.

    amx : The annotation matrix dataframe created by the load_amx function.

    elayers : The dataframe containing the 3D coordinates for the edges.

    color_col : The name of the column in the amx that contains the category
                that should be used for coloring the points in the plot.

    edge_index : The numerical index of the edge in the elayers dataframe
                 that is to be plotted.

    return : a Plotly figure object.
    """
    gene1 = elayers.iloc[edge_index]['Source']
    gene2 = elayers.iloc[edge_index]['Target']
    samples = elayers.iloc[edge_index]['Samples']

    # Generate the dataframe for the expression scatterplot
    sdata = pd.DataFrame(dict(x=gem.loc[gene1].values,y=gem.loc[gene2].values))
    sdata.index = gem.columns
    sdata = sdata.join(amx, how='left')

    # Generate the colors for the samples.
    categories = sdata[color_col].unique()
    num_categories = categories.shape[0]
    color_rep = pd.DataFrame({
        'Categories' : categories,
        'Color' : np.arange(0, num_categories)})
    sdata['color'] = sdata[color_col].replace(
                         to_replace=color_rep['Categories'].values,
                         value=color_rep['Color'].values)

    # Generate the Z-coordinates for the samples.
    color_rep = pd.DataFrame({
        'Categories' : categories,
        'Z' : np.arange(0, num_categories) / (num_categories - 1) - 0.5})
    sdata['z'] = sdata[color_col].replace(
                         to_replace=color_rep['Categories'].values,
                         value=color_rep['Z'].values)

    # Make samples in the cluster larger size
    sizes = pd.Series(list(samples))
    sizes = sizes.replace(to_replace=r'[^1]', value='5', regex=True)
    sizes = sizes.replace({'1': '10'})
    sizes = sizes.astype('int')
    fig2 = go.Figure(data=[go.Scatter3d(
                    x=sdata['x'],
                    y=sdata['y'],
                    z=sdata['z'],
                    mode='markers',
                    marker=dict(symbol='circle',size=sizes,
                                color=sdata['color'], colorscale='Jet'),
                    text= sdata['Sample'],
                    hoverinfo='text')])

    fig2.update_layout(
        height=600,
        title="3D Edge Co-Expression Scatterplot",
        showlegend=False,
        margin=dict(l=10, r=10, t=30, b=10),
        scene=dict(
          aspectmode="cube",
          xaxis=dict(showbackground=False, showline=True, zeroline=True, showgrid=True,
                     showticklabels=True, title=gene1,
                     showspikes=False),
          yaxis=dict(showbackground=False, showline=True, zeroline=True, showgrid=True,
                     showticklabels=True, title=gene2,
                     showspikes=False),
          zaxis=dict(showbackground=True, showline=True, zeroline=True, showgrid=True,
                     showticklabels=True, title='',
                     range=[-1.5, 1.5],
                     tickvals=[-1,1],
                     ticktext=['Sun', 'Shade'], showspikes=False),
        ),
        hovermode='closest',
        annotations=[dict(showarrow=False,
                        text="",
                        xref='paper', yref='paper',
                        x=0, y=0.1, xanchor='left', yanchor='bottom', font=dict(size=14))
                    ],
    )

    fig2.layout.scene.camera.projection.type = "orthographic"
    fig2.layout.scene.camera.eye = dict(x=0., y=0., z=-2)

    return fig2





def launch_application(net, gem, amx, vlayers, elayers, color_col):
    """
    Creates the Dash application.

    The Dash application will provide all of the interactive plots, tables and
    filters to interacitvely exploring the network.

    net :  The network dataframe created by the load_network function.

    gem :  The GEM dataframe created by the load_gem function.

    amx : The annotation matrix dataframe created by the load_amx function.

    vlayers : The dataframe containing the 3D coordinates for the nodes.

    elayers : The dataframe containing the 3D coordinates for the edges.

    """
    styles = {
        'pre': {
            'border': 'thin lightgrey solid',
            'overflowX': 'scroll'
        }
    }

    app = dash.Dash()
    app.layout = html.Div([
        html.Div(className='row', id = "header", children=[
            html.Img(
               src="https://raw.githubusercontent.com/SystemsGenetics/KINC/master/docs/images/kinc.png",
               style={
                  "height" : "55px",
                  "display" : "inline-block",
                  "padding" : "0px",
                  "margin" : "0px 10px 0px 10px"
               }),
            html.H1(
               children="3D Network Explorer",
               style={
                  "display" : "inline-block",
                  "padding" : "10px 0px 0px 0px",
                  "margin" : "0px",
                  "vertical-align" : "top"
               }),
            ]
        ),
        html.Div(className='row', children=[
            dcc.Graph(
              id = 'network-3dview',
              figure = create_network_plot(net, vlayers, elayers),
            )],
            style = {"width" : "45%",
                     "display" : "inline-block",
                     "border" : "1px solid black",
                     "padding" : "10px",
                     "background-color" : "black",
                     "margin" : "10px"},
        ),
        html.Div(className='row', children=[
            dcc.Graph(
              id = 'edge-expression-3dview',
              figure = create_expression_scatterplot(gem, amx, elayers, color_col),
            )],
            style = {"width" : "45%",
                     "display" : "inline-block",
                     "border" : "1px solid black",
                     "padding" : "10px",
                     "margin" : "10px"},
        ),
        # html.Div(className='row', children=[
        #   html.Div([
        #       dcc.Markdown("""
        #           **Click Data**
        #
        #           Click on points in the graph.
        #       """),
        #       html.Pre(id='click-data', style=styles['pre']),
        #   ], className='three columns'),
        # ])
    ])

    @app.callback(
        dash.dependencies.Output('edge-expression-3dview', 'figure'),
        [dash.dependencies.Input('network-3dview', 'clickData')])
    def update_expression_plot(clickData):
        if (clickData):
            points = clickData['points']
            found = re.match('^Edge: (.*?) \(co\) (.*?)$', points[0]['text'])
            if (found):
                edge_index = points[0]['pointNumber']
                figure = create_expression_scatterplot(gem, amx, elayers, color_col, edge_index)
                return(figure)

        figure = create_expression_scatterplot(gem, amx, elayers, color_col)
        return(figure)

    # @app.callback(
    #     dash.dependencies.Output('click-data', 'children'),
    #     [dash.dependencies.Input('network-3dview', 'clickData')])
    # def update_expression_plot(clickData):
    #     return json.dumps(clickData, indent=2)

    app.run_server(debug=True)





def main():
    """
    The main function.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--net', dest='net_path', type=str, required=True)
    parser.add_argument('--gem', dest='gem_path', type=str, required=True)
    parser.add_argument('--amx', dest='amx_path', type=str, required=True)
    parser.add_argument('--color_col', dest='color_col', type=str, required=True)
    args = parser.parse_args()

    # Load the input data.
    net = load_network(args.net_path)
    gem = load_gem(args.gem_path)
    amx = load_amx(args.amx_path)

    # Set the color column to use:
    color_col = args.color_col


    # Calculate a 2D layout for the network
    glayout = calculate_2d_layout(net)

    # Calculate the Z-coorinate positions for the verticies and edges.
    vlayers = get_vertex_zlayers(net, glayout)
    elayers = get_edge_zlayers(net, glayout)

    # Launch the dash application
    launch_application(net, gem, amx, vlayers, elayers, color_col)

    exit(0)



if __name__ == "__main__":
    main()
