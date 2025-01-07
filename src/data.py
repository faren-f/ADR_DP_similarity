import os
import pandas as pd
import numpy as np
import torch
from torch_geometric.data import HeteroData

class GraphConstruction ():
    def __init__(self, data_dir, graph_data_dir):

        # Read edge index
        edge_index_drug_adr = self._read_edge_index(data_dir, graph_data_dir['edge_index']['drug_adr'])
        edge_index_gene_drug = self._read_edge_index(data_dir, graph_data_dir['edge_index']['gene_drug'])
        edge_index_disease_dp = self._read_edge_index(data_dir, graph_data_dir['edge_index']['disease_dp'])                                        
        edge_index_gene_disease = self._read_edge_index(data_dir, graph_data_dir['edge_index']['gene_disease'])                                          
        edge_index_gene_gene = self._read_edge_index(data_dir, graph_data_dir['edge_index']['gene_gene'])
        
        # Read node attributes
        if os.path.exists(os.path.join(data_dir, graph_data_dir['node_attr']['drug']+'.csv')):
            node_attr_drug = self._read_node_attr(data_dir, graph_data_dir['node_attr']['drug'])
        else:
            Max_node_drug = torch.max(edge_index_drug_adr[0,:])
            node_attr_drug = torch.arange(0, Max_node_drug+1)


        if os.path.exists(os.path.join(data_dir, graph_data_dir['node_attr']['adr']+'.csv')):
            node_attr_adr = self._read_node_attr(data_dir, graph_data_dir['node_attr']['adr'])
        else:
            Max_node_adr = torch.max(edge_index_drug_adr[1,:])
            node_attr_adr = torch.arange(0, Max_node_adr+1)


        if os.path.exists(os.path.join(data_dir, graph_data_dir['node_attr']['disease']+'.csv')):
            node_attr_disease = self._read_node_attr(data_dir, graph_data_dir['node_attr']['disease'])
        else:
            Max_node_disease = torch.max(edge_index_disease_dp[0,:])
            node_attr_disease = torch.arange(0, Max_node_disease+1)


        if os.path.exists(os.path.join(data_dir, graph_data_dir['node_attr']['dp']+'.csv')):
            node_attr_dp = self._read_node_attr(data_dir, graph_data_dir['node_attr']['dp'])
        else:
            Max_node_dp = torch.max(edge_index_disease_dp[1,:])
            node_attr_dp = torch.arange(0, Max_node_dp+1)


        if os.path.exists(os.path.join(data_dir, graph_data_dir['node_attr']['gene']+'.csv')):
            node_attr_gene = self._read_node_attr(data_dir, graph_data_dir['node_attr']['gene'])
        else:
            Max_node_gene = torch.max(torch.max(edge_index_gene_gene[0,:]),  torch.max(edge_index_gene_gene[1,:]))
            node_attr_gene = torch.arange(0, Max_node_gene+1)



        self.num_nodes_dict = {'adr': torch.tensor(node_attr_adr.shape[0]), 'dp': torch.tensor(node_attr_dp.shape[0]),
                          'drug': torch.tensor(node_attr_drug.shape[0]), 'disease': torch.tensor(node_attr_disease.shape[0]),
                          'gene': torch.tensor(node_attr_gene.shape[0])}



        self.edge_index_drug_adr = edge_index_drug_adr
        self.edge_index_gene_drug = edge_index_gene_drug
        self.edge_index_disease_dp = edge_index_disease_dp                                        
        self.edge_index_gene_disease = edge_index_gene_disease                                          
        self.edge_index_gene_gene = edge_index_gene_gene

        self.node_attr_drug = node_attr_drug
        self.node_attr_adr = node_attr_adr
        self.node_attr_disease = node_attr_disease
        self.node_attr_dp = node_attr_dp
        self.node_attr_gene = node_attr_gene
                          
    def _read_edge_index(self, data_dir, file_name):
        edge_index = pd.read_csv(os.path.join(data_dir, file_name+'.csv'))
        edge_index = np.array(edge_index)
        edge_index = edge_index[:, [2,3]]
        edge_index = edge_index.astype(np.int64)
        edge_index = np.transpose(edge_index)
        edge_index = torch.tensor(edge_index, dtype = torch.long)
        return edge_index
    
    def _read_node_attr (self, data_dir, file_name):
        node_attr = pd.read_csv(os.path.join(data_dir, file_name+'.csv'))
        node_attr = np.array(node_attr)
        node_attr = torch.tensor(node_attr, dtype = torch.float)
        return node_attr

    def create_graph(self):
        data = HeteroData()
        data.num_nodes_dict = self.num_nodes_dict
        data['adr'].x = self.node_attr_adr
        data['dp'].x = self.node_attr_dp
        data['drug'].x = self.node_attr_drug
        data['disease'].x = self.node_attr_disease
        data['gene'].x = self.node_attr_gene
        
        data['drug', 'conntectedTo', 'adr'].edge_index = self.edge_index_drug_adr
        data['disease', 'conntectedTo', 'dp'].edge_index = self.edge_index_disease_dp
        data['gene', 'targetedBy', 'drug'].edge_index = self.edge_index_gene_drug
        data['gene', 'associatedWith', 'disease'].edge_index = self.edge_index_gene_disease
        data['gene', 'interactWith', 'gene'].edge_index = self.edge_index_gene_gene
        return data


