import os
import torch
from torch.nn import Embedding
import torch.nn.functional as F 
from torch.nn import Linear

import torch_geometric
from torch_geometric.nn import to_hetero, GraphConv
from torch_geometric.nn import MLP

class GNNEncoder(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GraphConv((-1,-1), hidden_channels)
        self.conv2 = GraphConv((-1,-1), hidden_channels)
        #self.conv3 = GraphConv((-1,-1), hidden_channels)
        self.conv4 = GraphConv((-1,-1), out_channels)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()
        #x = self.conv3(x, edge_index).relu()
        x = self.conv4(x, edge_index)
        return(x)

class EdgeDecoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.linear = Linear(1,1, bias=True)

    def forward(self, z_dict, edge_label_index, type, apply_linear=False):
        row_0, row_1 = edge_label_index
        

        if type == 1:
            sim = F.cosine_similarity(z_dict['drug'][row_0], z_dict['adr'][row_1], dim = -1, eps = 1e-6)
            # sim = (z_dict['adr'][row_0] * z_dict['drug'][row_1]).sum(dim=-1)  # dot product 

        elif type == 2:
            sim = F.cosine_similarity(z_dict['gene'][row_0], z_dict['drug'][row_1], dim = -1, eps = 1e-6)
            # sim = (z_dict['drug'][row_0] * z_dict['gene'][row_1]).sum(dim=-1)  # dot product 


        elif type == 3:
            sim = F.cosine_similarity(z_dict['disease'][row_0], z_dict['dp'][row_1], dim = -1, eps = 1e-6)
            # sim = (z_dict['dp'][row_0] * z_dict['disease'][row_1]).sum(dim=-1)  # dot product 

        elif type == 4:
            sim = F.cosine_similarity(z_dict['gene'][row_0], z_dict['disease'][row_1], dim = -1, eps = 1e-6)
            # sim = (z_dict['disease'][row_0] * z_dict['gene'][row_1]).sum(dim=-1)  # dot product

        elif type == 5:
            sim = F.cosine_similarity(z_dict['gene'][row_0], z_dict['gene'][row_1], dim = -1, eps = 1e-6)

        else:
            print('type must be <= 5')
            
        if apply_linear:
            y = self.linear(sim.view(-1, 1))
        return y.squeeze() if apply_linear else sim


class Model(torch.nn.Module):
    def __init__(self, num_nodes_dict, hidden_channels, data):
        super().__init__()
        self.data = data
        # Initialize a separate embedding layer for each node type
        self.embeddings = torch.nn.ModuleDict({
            node_type: Embedding(num_nodes, hidden_channels) for node_type, num_nodes in num_nodes_dict.items()
        })
        
        self.encoder = GNNEncoder(hidden_channels, hidden_channels)
        self.encoder = to_hetero(self.encoder, data.metadata(), aggr='sum')
        self.decoder = EdgeDecoder()

    def forward (self, x_dict, edge_index_dict, edge_label_index, type):
        
        emb_adr = self.embeddings['adr'](x_dict['adr'])
        emb_dp = self.embeddings['dp'](x_dict['dp'])
        emb_drug = self.embeddings['drug'](x_dict['drug'])
        emb_disease = self.embeddings['disease'](x_dict['disease'])
        emb_gene = self.embeddings['gene'](x_dict['gene'])

        emb_dict = {
            'adr': emb_adr, 
            'dp': emb_dp,
            'drug': emb_drug, 
            'disease': emb_disease,
            'gene': emb_gene  
        }
        
        z_dict = self.encoder(emb_dict, edge_index_dict)
       
        return self.decoder(z_dict, edge_label_index, type)

    
    def forward_encoder(self, x_dict, edge_index_dict, raw_data_dir, graph_data_dir):
        
        # if we have node features for a node or not
        if os.path.exists(os.path.join(raw_data_dir, graph_data_dir['node_attr']['adr']+'.csv')):
            emb_adr = self.data['adr'].x
        else:
            emb_adr = self.embeddings['adr'](x_dict['adr'])

        if os.path.exists(os.path.join(raw_data_dir, graph_data_dir['node_attr']['dp']+'.csv')):
            emb_dp = self.data['dp'].x
        else:
            emb_dp = self.embeddings['dp'](x_dict['dp'])

        if os.path.exists(os.path.join(raw_data_dir, graph_data_dir['node_attr']['drug']+'.csv')):
            emb_drug = self.data['drug'].x
        else:
            emb_drug = self.embeddings['drug'](x_dict['drug'])


        if os.path.exists(os.path.join(raw_data_dir, graph_data_dir['node_attr']['disease']+'.csv')):
            emb_disease = self.data['disease'].x
        else:
            emb_disease = self.embeddings['disease'](x_dict['disease'])


        if os.path.exists(os.path.join(raw_data_dir, graph_data_dir['node_attr']['gene']+'.csv')):
            emb_gene = self.data['gene'].x
        else:
            emb_gene = self.embeddings['gene'](x_dict['gene'])
            

        emb_dict = {
            'adr': emb_adr, 
            'dp': emb_dp,
            'drug': emb_drug, 
            'disease': emb_disease,
            'gene': emb_gene  
        }
        z_dict = self.encoder(emb_dict, edge_index_dict)
        return z_dict




