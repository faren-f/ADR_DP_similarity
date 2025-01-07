import pandas as pd
import torch
from sklearn.metrics import roc_auc_score



class Evaluation():
    def __init__(self, y_true, y_pred):
        self.y_true = y_true
        self.y_pred = y_pred
        self.y_pred_pos = y_pred[y_true == 1]
        self.y_pred_neg = y_pred[y_true == 0]

    def eval(self, k_list):

        df = pd.DataFrame(columns = ['prc@' + str(k) for k in k_list] + 
                                    ['hits@' + str(k) for k in k_list] +
                                    ['mrr_hits@' + str(k) for k in [1,3,10,20,50,100]] +
                                    ['mrr', 'mrr2', 'auc'])
        Precision_k_list = []
        Hits_k_list = []
        for k in k_list:
            # Precision_k
            Precision_k_list.append(self.precision_k(k))
            Hits_k_list.append(self.hits_k(k))

        MRR_list = self.mrr()
        MRR2 = [self.mrr2()]
        AUC = [self.auc()]
        df.loc[0] = Precision_k_list + Hits_k_list + MRR_list + MRR2 + AUC
        return df

    def precision_k(self, k):
        """
        the fraction of true links that appear in the first 𝑘 link of the sorted rank list.
        
        y_true: A tensor of ground truth (0 and 1).
        y_pred: A tensor of logits.
        k: Number of top elements to look at for computing precision.
        """
        # Get indices of the predictions with the highest scores
        if k > len(self.y_pred):
            k = len(self.y_pred)
        
        topk_indices = torch.topk(self.y_pred, k).indices
        
        # Calculate precision
        value = self.y_true[topk_indices].float().mean()
        return value.item()
    
    
    def hits_k(self, k):
        
        if k > len(self.y_pred_neg):
            return(1)
            
        # find the kth score in the negative predictions as the threshold
        kth_score_in_negative_edges = torch.topk(self.y_pred_neg, k=k)[0][-1]
        # Count the number of positive predictions > the threshold the / len(y_pred_pos)
        hitsK = float(torch.sum(self.y_pred_pos > kth_score_in_negative_edges).cpu()) / len(self.y_pred_pos)
        return hitsK
    
    
    def mrr(self):
        '''
            compute mrr
            y_pred_neg is an array with shape (batch size, num_entities_neg).
            y_pred_pos is an array with shape (batch size, num_entities_pos)
        '''
        
        # calculate ranks
        y_pred_pos = self.y_pred_pos.view(-1, 1)
        # optimistic rank: "how many negatives have at least the positive score?"
        # ~> the positive is ranked first among those with equal score
        optimistic_rank = (self.y_pred_neg >= y_pred_pos).sum(dim=1)
        # pessimistic rank: "how many negatives have a larger score than the positive?"
        # ~> the positive is ranked last among those with equal score
        pessimistic_rank = (self.y_pred_neg > y_pred_pos).sum(dim=1)
        ranking_list = 0.5 * (optimistic_rank + pessimistic_rank) + 1
    
        hits1_list = (ranking_list <= 1).to(torch.float)
        hits3_list = (ranking_list <= 3).to(torch.float)
        hits10_list = (ranking_list <= 10).to(torch.float)
        hits20_list = (ranking_list <= 20).to(torch.float)
        hits50_list = (ranking_list <= 50).to(torch.float)
        hits100_list = (ranking_list <= 100).to(torch.float)
        mrr_list = 1./ranking_list.to(torch.float)
    
        return [hits1_list.mean().item(),
                hits3_list.mean().item(),
                hits10_list.mean().item(),
                hits20_list.mean().item(),
                hits50_list.mean().item(),
                hits100_list.mean().item(),
                mrr_list.mean().item()]
    
    def mrr2(self):
        """
        Calculate the Mean Reciprocal Rank (MRR) for link prediction.
    
        y_true: A tensor of ground truth (0 and 1).
        y_pred: A tensor of logits.
        """
    
        # Sort y_pred in descending order and get the indices
        sorted_indices = torch.argsort(self.y_pred, descending=True) +1
    
        # Calculate the mean reciprocal rank
        mrr = (1 / sorted_indices[self.y_true == 1]).mean()
        return mrr.item()
    
    
    def auc(self):
        # try:
        auc_score = roc_auc_score(self.y_true.detach().cpu().numpy(), self.y_pred.detach().cpu().numpy())
        # except ValueError:  # Handle case with only one class present in y_true
        #    auc_score = float('nan')
        return auc_score





