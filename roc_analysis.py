import sys
import matplotlib.pyplot as plt
from sklearn import metrics

class DrawPlot:
    def __init__(self):
        plt.figure()
        lw = 2
        files = ['./tpr_fpr_hotnet2.txt', './tpr_fpr_our.txt']
        #colors = ['aqua', 'darkorange', 'cornflowerblue', 'red', 'purple', 'brown', 'green', 'navy', 'yellow', 'grey']
        colors = ['aqua', 'purple']
        labels = ['H', 'O']
        for file,color,lab in zip(files,colors,labels):
            with open(file) as cur_file:
                tpr = []
                fpr = []
                for line in cur_file:
                    cols = line.split()
                    tpr.append(float(cols[0]))
                    fpr.append(float(cols[1]))
            plt.plot(fpr, tpr, color=color, lw=lw, label='%s (%f )' % (lab,metrics.auc(fpr,tpr)))
            plt.plot(fpr, tpr, color=color, lw=lw, label='%s' % lab)

        plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curves of Employed Measures')
        plt.legend(loc="lower right")
        plt.show()

if __name__ == "__main__":
    inputdata = DrawPlot()
