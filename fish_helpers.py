import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from metadata import Metadata
import os
from fish_results import HybeData
from collections import defaultdict, Counter
from scipy import stats
import seaborn as sns
from scipy.stats import spearmanr

def ReadsPerGene_fun(system='cornea'):
    if system=='cornea':
        f = '/bigstore/GeneralStorage/Zach/Cornea_RNAseq/Aligned/ReadsPerGene.xlsx'
        ReadsPerGene = pd.read_excel(f)
        ReadsPerGene.index = ReadsPerGene.GeneIDs
    elif system=='3t3':
        f = '/bigstore/GeneralStorage/Evan/NFKB_MERFISH/Calibration_Set_Data/3T3_Calibration_Set/RNA_Seq_Data/Hoffmann_IFNAR_KO_3T3_TNF.txt'
        ReadsPerGene = pd.read_csv(f,sep='\t')
        gids = []
        for gid in ReadsPerGene.Geneid:
            gids.append(gid.split('.')[0])
        ReadsPerGene.Geneid = gids
        ReadsPerGene.index = gids
    else:
        print('Unknown System')
    return ReadsPerGene

def merfish_correlation(spotcalls,color='b',alpha=1,label='spotcalls',system='cornea',ReadsPerGene=False):
    import pandas as pd
    if not isinstance(ReadsPerGene,pd.core.frame.DataFrame):
        if system=='cornea':
            f = '/bigstore/GeneralStorage/Zach/Cornea_RNAseq/Aligned/ReadsPerGene.xlsx'
            ReadsPerGene = pd.read_excel(f)
            ReadsPerGene.index = ReadsPerGene.GeneIDs
        elif system=='3t3':
            f = '/bigstore/GeneralStorage/Evan/NFKB_MERFISH/Calibration_Set_Data/3T3_Calibration_Set/RNA_Seq_Data/Hoffmann_IFNAR_KO_3T3_TNF.txt'
            ReadsPerGene = pd.read_csv(f,sep='\t')
            gids = []
            for gid in ReadsPerGene.Geneid:
                gids.append(gid.split('.')[0])
            ReadsPerGene.Geneid = gids
            ReadsPerGene.index = gids
        else:
            print('Unknown System')
    GeneList = pd.read_excel('/bigstore/GeneralStorage/Zach/MERFISH/Inflammatory/InflammationGeneList.xlsx')
    GeneList.index = GeneList.Gene
    counts = []
    fpkms = []
    from collections import defaultdict, Counter
    FISH_Spots = Counter(spotcalls.gene)
    for gn,cc in FISH_Spots.items():
        if 'blank' in gn:
            continue
        else:
            gid = GeneList.loc[gn]['Gene_ID']
            if system=='cornea':
                reads = ReadsPerGene.loc[gid].Unstranded
            elif system=='3t3':
                reads = ReadsPerGene.loc[gid]['IFNAR-TNF-3_tot.bam']
            else:
                print('Unknown System')
            fpkm = reads/GeneList.loc[gn]['Length']
            if isinstance(fpkm,np.float64):
                if cc<2:
                    continue
                counts.append(cc)
                fpkms.append(fpkm)
    from scipy.stats import spearmanr
    import matplotlib.pyplot as plt
    plt.scatter(np.log10(fpkms),np.log10(counts),c=color,alpha=alpha,label=label)
    print(spearmanr(fpkms,counts))
    plt.suptitle('Untrimmed FPKM vs Spot Count')
    plt.ylabel('log10 MERFISH Spot Count')
    plt.xlabel('log10 RNAseq FPKM')
    plt.legend()

def ptl_hist(spotcalls,bins=1000,colors='rkb',alpha=0.5,column='ave'):
    import matplotlib.pyplot as plt
    plt.hist(np.log10(spotcalls[column]),bins=bins,color=colors[0],alpha=alpha,label='raw')
    plt.hist(np.log10(spotcalls[spotcalls.npixels>1][column]),bins=bins,color=colors[1],alpha=alpha,label='npixels>1')
    plt.hist(np.log10(spotcalls[spotcalls.ssum>2**12][column]),bins=bins,color=colors[2],alpha=alpha,label='ssum>2**12')
    plt.ylabel('Counts')
    plt.xlabel(str('Log10 '+str(column)))
    plt.legend()
    plt.show()
    
def spotcalls_qc(spotcalls,system='cornea',colors='rkb',alpha=0.5,bins=1000):
    import matplotlib.pyplot as plt
    ReadsPerGene=ReadsPerGene_fun(system='cornea')
    merfish_correlation(spotcalls,color=colors[0],alpha=alpha,label='raw',system='cornea',ReadsPerGene=ReadsPerGene)
    merfish_correlation(spotcalls[spotcalls.npixels>1],color=colors[1],alpha=alpha,label='npixels>1',system='cornea',ReadsPerGene=ReadsPerGene)
    merfish_correlation(spotcalls[spotcalls.ssum>2**12],color=colors[2],alpha=alpha,label='ssum>2**12',system='cornea',ReadsPerGene=ReadsPerGene)
    plt.show()
    ptl_hist(spotcalls,bins=bins,colors=colors,alpha=alpha,column='ave')
    ptl_hist(spotcalls,bins=bins,colors=colors,alpha=alpha,column='ssum')
    
def onfly_qc(md,path=False):
    if path==True:
        from metadata import Metadata
        md = Metadata(md)
    import pickle
    import os
    import time
    i=0
    while i==0:
        try:
            tforms = pickle.load(open(os.path.join(md_path,'results','tforms.pkl'),'rb'))
            beads = pickle.load(open(os.path.join(md_path,'results','beads.pkl'),'rb'))
            i=1
        except:
            time.sleep(5)
    print(len(tforms['good'].keys()),' Good Positions')
    print(len(tforms['bad'].keys()),' Failed Positions')
    import matplotlib.pyplot as plt
    X = []
    Y = []
    good_pos = tforms['good'].keys()
    bad_pos = tforms['bad'].keys()
    for pos in good_pos:
        x,y = md.image_table[md.image_table.Position==pos].XY.iloc[0]
        X.append(x)
        Y.append(y)
    plt.scatter(X,Y,c='g',label='good')
    X = []
    Y = []
    for pos in bad_pos:
        if pos in good_pos:
            continue
        x,y = md.image_table[md.image_table.Position==pos].XY.iloc[0]
        X.append(x)
        Y.append(y)
    plt.scatter(X,Y,c='r',label='bad')
    plt.xlabel('X Stage Coordinate')
    plt.ylabel('Y Stage Coordinate')
    plt.legend()
    plt.show()
    
def photobleach_qc(md,path=True,pos=False):
    import matplotlib.pyplot as plt
    if path==True:
        from metadata import Metadata
        md = Metadata(md)
    if pos ==False:
        pos = md.image_table.Position.iloc[0]
    for acq in md.image_table[md.image_table.Position==pos].acq.unique():
        if 'hybe' in acq:
            stk = md.stkread(Position=pos,Channel='FarRed',acq=acq)
            plt.plot(range(stk.shape[2]),np.mean(np.mean(stk,axis=0),axis=0),label=acq)
    plt.title('FarRed')
    plt.xlabel('Z index')
    plt.ylabel('Average Intensity')
    plt.legend()
    plt.show()
    for acq in md.image_table[md.image_table.Position==pos].acq.unique():
        if 'hybe' in acq:
            stk = md.stkread(Position=pos,Channel='Orange',acq=acq)
            plt.plot(range(stk.shape[2]),np.mean(np.mean(stk,axis=0),axis=0),label=acq)
    plt.title('Orange')
    plt.xlabel('Z index')
    plt.ylabel('Average Intensity')
    plt.legend()
    plt.show()

def img2stage_coordinates(spotcalls,md,path=False,pixelsize=0.109,cameradirection=[1,1],verbose=False):
    if path==True:
            from metadata import Metadata
            md = Metadata(md)
    X = []
    Y = []
    for pos in spotcalls.posname.unique():
        if verbose==True:
            print(pos)
        coordX,coordY = md.image_table[md.image_table.Position==pos].XY.iloc[0]
        pos_temp = spotcalls[spotcalls.posname==pos]
        pos_centroid = pos_temp.centroid
        for xy in pos_centroid:
            rx = xy[1]-2048
            ry = xy[0]#-2048
            x = rx*pixelsize*cameradirection[0]+coordX
            y = ry*pixelsize*cameradirection[1]+coordY
            X.append(x)
            Y.append(y)
    spotcalls['CoordX'] = X
    spotcalls['CoordY'] = Y
    return spotcalls

def progress_update(md,path=False):
    if path==True:
        from metadata import Metadata
        md = Metadata(md)
    finished = []
    non_finished = []
    import os
    for pos in os.listdir(os.path.join(md.base_pth,'codestacks')):
        try:
            processed = pickle.load(open(os.path.join(md.base_pth,'codestacks',pos,'processing.pkl'),'rb'))
        except:
            processed=[]
        if len(processed)==18:
            finished.append(pos)
        else:
            non_finished.append(pos)
    print(len(finished),' Positions Finished')
    print(len(non_finished),' Positions Not Finished')