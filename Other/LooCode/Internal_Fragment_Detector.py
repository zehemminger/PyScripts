#!/usr/bin/env python

import os
import sys
import csv
import time
import math
import pickle
import argparse
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import multiprocessing
from array import array
from functools import partial
import matplotlib.pyplot as plt
from matplotlib import transforms


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence_path", type=str, help="Path to txt file containing protein sequence.")
    parser.add_argument("fragmentation_path", type=str, help="Path to excel file with detected peaks")
    parser.add_argument("-b", "--base_path", type=str, dest="base_path",default='', help="Path to folder containing AminoAcids.csv and Modifications.csv (default is current location)")
    parser.add_argument("-o", "--out_path", type=str, dest="out_path",default='', help="Path to folder where you want the results saved (default is current location)")
    parser.add_argument("-p", "--nthreads", type=int, dest="ncpu", default=0, action='store', help="Number of cores to utilize (default uses all threads available).")
    args = parser.parse_args()
    
def generate_column(frag_type,out_path):
    seq_df = pd.read_csv(os.path.join(out_path,'seq_df.csv'),index_col=0)
    frag_list = []
    seq_list = []
    total_length = len(seq_df.index)
    if frag_type =='N':
        for length in seq_df.index:
            frag_list.append(np.sum(seq_df['MSs'].loc[0:length]))
            seq_list.append(''.join(i for i in seq_df['AAs'].loc[0:length]))
    elif frag_type =='C':
        for length in seq_df.index:
            frag_list.append(np.sum(seq_df['MSs'].loc[total_length-length:]))
            seq_list.append(''.join(i for i in seq_df['AAs'].loc[total_length-length:]))
        frag_list.reverse()
        seq_list.reverse()
    elif 'I' in frag_type:
        length = int(frag_type.split('_')[1])
        for start_site in seq_df.index:
            if start_site==0:
                frag_list.append(np.nan) #Same as N term
                seq_list.append(np.nan)
            elif start_site==total_length:
                frag_list.append(np.nan) #Same as C term
                seq_list.append(np.nan)
            elif start_site+length>=total_length:
                frag_list.append(np.nan) #Doesnt Exist
                seq_list.append(np.nan)
            else:
                frag_list.append(np.sum(seq_df['MSs'].loc[start_site:start_site+length]))
                seq_list.append(''.join(i for i in seq_df['AAs'].loc[start_site:start_site+length]))
    return frag_type,frag_list,seq_list

def Lookup_Wrapper(Sequence,out_path,ncpu=4):
    seq_df = pd.read_csv(os.path.join(out_path,'seq_df.csv'),index_col=0)
    Frag_Result = seq_df.copy()
    Seq_Result = seq_df.copy()
    Input = ['N','C']
    for length in seq_df.index:
        Input.append('I_'+str(length))
    if 'linux' in sys.platform:
        with multiprocessing.Pool(ncpu) as p:
            sys.stdout.flush()
            pfunc = partial(generate_column,out_path=out_path)
            for frag_type,frag_list,seq_list in p.imap(pfunc,Input, chunksize=1):
                Frag_Result[frag_type] = frag_list
                Seq_Result[frag_type] = seq_list
            p.close()
            sys.stdout.flush()
    elif 'win' in sys.platform:
        if __name__ ==  '__main__':
            with multiprocessing.Pool(ncpu) as p:
                sys.stdout.flush()
                pfunc = partial(generate_column,out_path=out_path)
                for frag_type,frag_list,seq_list in p.imap(pfunc,Input, chunksize=1):
                    Frag_Result[frag_type] = frag_list
                    Seq_Result[frag_type] = seq_list
                p.close()
                sys.stdout.flush()
    return Frag_Result,Seq_Result

def match_peaks(Fragment,shift,mod,frag_array,seq_array):
    Fragment -= shift
    diff_array = np.abs(np.subtract(Fragment,frag_array))
    Error = np.nanmin(diff_array)
    x,y = np.where(diff_array==Error)
    seqs = [i+'_'+mod for i in seq_array[x,y]]
    return Error,seqs,x,y

def seq_mod1_mod2_to_peak(result,N_Frag_Mod,C_Frag_Mod,AminoAcids):
    seq,mod1,mod2 = result.split('_')
    peak = []
    for aa in seq:
        peak.append(AminoAcids.loc[aa])
    peak.append(N_Frag_Mod[mod1])
    peak.append(C_Frag_Mod[mod2])
    return np.sum(peak)

def calc_end_length(start_sites,Assignment):
    end_sites = []
    length = []
    for i,ss in enumerate(start_sites):
        l = len(Assignment[i].split('_')[0])
        length.append(l)
        end_sites.append(ss+l)
    return end_sites,length

def find_best_match(Fragment,out_path):
    N_Frag_Mod = pickle.load(open(os.path.join(out_path,'N_Frag_Mod.pkl'),'rb'))
    C_Frag_Mod = pickle.load(open(os.path.join(out_path,'C_Frag_Mod.pkl'),'rb'))
    frag_array = pickle.load(open(os.path.join(out_path,'frag_array.pkl'),'rb'))
    seq_array = pickle.load(open(os.path.join(out_path,'seq_array.pkl'),'rb'))
    AminoAcids = pd.read_csv(os.path.join(out_path,'AminoAcids.csv'),index_col=0)
    Errors = []
    Seqs = []
    X = []
    Y = []
    for mod1,shift1 in N_Frag_Mod.items():
        for mod2,shift2 in C_Frag_Mod.items():
            mod = mod1+'_'+mod2
            shift = shift1+shift2
            Error,seqs,x,y = match_peaks(Fragment,shift,mod,frag_array,seq_array)
            X.append(x)
            Y.append(y)
            Errors.append(Error)
            Seqs.append(seqs)
    min_Error = np.min(Errors)
    Assignment = Seqs[np.where(Errors==min_Error)[0][0]]
    calc_peaks = seq_mod1_mod2_to_peak(Assignment[0],N_Frag_Mod,C_Frag_Mod,AminoAcids)[0]
    start_sites = X[np.where(Errors==min_Error)[0][0]]
    end_sites,length = calc_end_length(start_sites,Assignment)
    return Fragment,min_Error,Assignment,calc_peaks,start_sites,end_sites,length

def match_wrapper(Fragments,seq_df,Sequence_df,out_path,ncpu=4):
    Fragment_list = []
    min_Errors = []
    Assignments = []
    Calculated = []
    Start_sites = []
    End_sites = []
    Length = []
    N_Type = []
    C_Type = []
    Input = Fragments['Fragment']
    with multiprocessing.Pool(ncpu) as p:
        sys.stdout.flush()
        pfunc = partial(find_best_match,out_path=out_path)
        for Fragment,min_Error,Assignment,calc_peaks,start_sites,end_sites,length in p.imap(pfunc,Input, chunksize=1):
            Fragment_list.append(Fragment)
            min_Errors.append(min_Error)
            Assignments.append(Assignment)
            Calculated.append(calc_peaks)
            Start_sites.append(start_sites)
            End_sites.append(end_sites)
            Length.append(length)
            N_Type.append([i.split('_')[1] for i in Assignment])
            C_Type.append([i.split('_')[2] for i in Assignment])
        p.close()
        sys.stdout.flush()
    Results = pd.DataFrame()
    Results['Fragment'] = Fragment_list
    Results['Error'] = min_Errors
    Results['Assignment'] = Assignments
    Results['Calculated'] = Calculated
    Results['Start_site'] = Start_sites
    Results['End_site'] = End_sites
    Results['Length'] = Length
    Results['N_Type'] = N_Type
    Results['C_Type'] = C_Type
    return Results

def rainbow_text(x, y, strings, colors, orientation='horizontal',
                 ax=None, **kwargs):
    """
    Take a list of *strings* and *colors* and place them next to each
    other, with text strings[i] being shown in colors[i].

    Parameters
    ----------
    x, y : float
        Text position in data coordinates.
    strings : list of str
        The strings to draw.
    colors : list of color
        The colors to use.
    orientation : {'horizontal', 'vertical'}
    ax : Axes, optional
        The Axes to draw into. If None, the current axes will be used.
    **kwargs
        All other keyword arguments are passed to plt.text(), so you can
        set the font size, family, etc.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    assert orientation in ['horizontal', 'vertical']
    if orientation == 'vertical':
        kwargs.update(rotation=90, verticalalignment='bottom')

    for s, c in zip(strings, colors):
        text = ax.text(x, y, s + " ", color=c, transform=t, **kwargs)

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        if orientation == 'horizontal':
            t = transforms.offset_copy(
                text.get_transform(), x=ex.width, units='dots')
        else:
            t = transforms.offset_copy(
                text.get_transform(), y=ex.height, units='dots')
def view_coverage(Results,Sequence,Colors = ['purple','blue','red'],figsize=[6,6],font_size=18,step=25,shift=1):
    Frag_df = pd.DataFrame()
    Frag_df['Start_Frag_sites'] = list(itertools.chain.from_iterable(Results['Start_site']))
    Frag_df['End_Frag_sites'] = list(itertools.chain.from_iterable(Results['End_site']))
    Total_frag_loc = np.unique(Frag_df)
    print(len(Total_frag_loc),' Total Fragments')
    Terminal_frag_loc = np.unique(Frag_df[(Frag_df['Start_Frag_sites']==1) | (Frag_df['End_Frag_sites']==len(Sequence))])
    print(len(Terminal_frag_loc),' Terminal Fragments')
    Internal_frag_loc = np.unique(Frag_df[(Frag_df['Start_Frag_sites']!=1) & (Frag_df['End_Frag_sites']!=len(Sequence))])
    print(len(Internal_frag_loc),' Internal Fragments')
    Shared_frag_loc = np.unique([i for i in Terminal_frag_loc if i in Internal_frag_loc])
    print(len(Shared_frag_loc),' Shared Fragments')
    plt.figure(figsize=figsize)
    length = len(Sequence)+len(Total_frag_loc)
    step = step
    lines = math.ceil(len(Sequence)/step)
    temp_words = []
    temp_colors = []
    shift = shift
    line = 0
    size = font_size
    for i,aa in enumerate(Sequence):
        temp_words.append(aa)
        temp_colors.append('black')
        if i==len(Sequence)-shift:
            continue
        if len(Colors)==1:
            if i+shift in Total_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[0]))
            else:
                temp_words.append('|')
                temp_colors.append('white')
        elif len(Colors)==2:
            if i+shift in Shared_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[0]))
            elif i+shift in Terminal_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[1]))
            elif i+shift in Internal_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[0]))
            else:
                temp_words.append('|')
                temp_colors.append('white')
        elif len(Colors)==3:
            if i+shift in Shared_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[0]))
            elif i+shift in Terminal_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[1]))
            elif i+shift in Internal_frag_loc:
                temp_words.append('|')
                temp_colors.append(str(Colors[2]))
            else:
                temp_words.append('|')
                temp_colors.append('white')
        if (i+1)%step==0:
            rainbow_text(0, ((lines-line)/lines), temp_words, temp_colors, size=size)
            temp_words = []
            temp_colors = []
            line+=1
    rainbow_text(0, ((lines-line)/lines), temp_words, temp_colors, size=size)
    plt.axis('off')
    plt.show()
    
if __name__ == '__main__':
    Master_Start = time.time()
#     global seq_df,AminoAcids
#     global N_Frag_Mod,C_Frag_Mod,frag_array,seq_array
    
    #Atom Masses
    n = 14.003074
    o = 15.994915
    h = 1.007825
    p = 1.007276
    car = 12.000000

    sequence_path = args.sequence_path
    out_path = args.out_path
    base_path = args.base_path
    fragmentation_path = args.fragmentation_path
    ncpu = args.ncpu
    print(args)
    
    if base_path == '':
        base = os.getcwd()
    else:
        base = base_path
    if out_path == '':
        out_path = os.getcwd()
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    if ncpu==0:
        ncpu = multiprocessing.cpu_count()
    #Import Amino Acid Sequence
    with open(sequence_path,"r") as Seq:
        Sequence = Seq.read()
    AminoAcids = pd.read_csv(os.path.join(base,'AminoAcids.csv'),header=None,index_col=0,usecols=[0,1],names=['AAs','MSs'])
    AminoAcids.to_csv(os.path.join(out_path,'AminoAcids.csv'))
    seq_AAs = []
    seq_MSs = []
    for aa in Sequence:
        seq_AAs.append(aa)
        seq_MSs.append(AminoAcids.loc[aa][0])
    seq_df = pd.DataFrame()
    seq_df['AAs'] = seq_AAs
    seq_df['MSs'] = seq_MSs
    seq_df.to_csv(os.path.join(out_path,'seq_df.csv'))
    Total_Mass_Two = np.sum(seq_df['MSs']) + o + 2 * h
    
    start = time.time()
    Frag_Result,Seq_Result = Lookup_Wrapper(Sequence,out_path,ncpu=ncpu)
    print(round(time.time()-start,2),'Seconds to generate lookup table')
    Seq_Result.to_csv(os.path.join(out_path,'Seq_Result.csv'))
    Frag_Result.to_csv(os.path.join(out_path,'Frag_Result.csv'))
    
    N_Frag_Mod = {}
    N_Frag_Mod[''] = 0
    N_Frag_Mod['a'] = - car - o
    N_Frag_Mod['ap'] = N_Frag_Mod['a']+p
    N_Frag_Mod['b'] = 0
    N_Frag_Mod['bp'] = N_Frag_Mod['b']+p
    N_Frag_Mod['c'] = n + h + h + h
    N_Frag_Mod['cp'] = N_Frag_Mod['c']+p
    pickle.dump(N_Frag_Mod,open(os.path.join(out_path,'N_Frag_Mod.pkl'),'wb'))
    C_Frag_Mod = {}
    C_Frag_Mod[''] = 0
    C_Frag_Mod['x'] = o + car + o
    C_Frag_Mod['xp'] = C_Frag_Mod['x']+p
    C_Frag_Mod['y'] = o + h
    C_Frag_Mod['yp'] = C_Frag_Mod['y']+p+h
    C_Frag_Mod['z'] = o + h - n - h
    C_Frag_Mod['zp'] = C_Frag_Mod['z']+p+h
    pickle.dump(C_Frag_Mod,open(os.path.join(out_path,'C_Frag_Mod.pkl'),'wb'))         
    Fragments = pd.read_excel(fragmentation_path,header=None)
    Fragments.columns = ['Fragment']
    #Modifications = pd.read_csv(os.path.join(base,"Modifications.csv"),header=None)
    
    frag_array = np.array(Frag_Result.drop(columns='AAs'))
    seq_array = np.array(Seq_Result.drop(columns='MSs'))
    pickle.dump(frag_array,open(os.path.join(out_path,'frag_array.pkl'),'wb')) 
    pickle.dump(seq_array,open(os.path.join(out_path,'seq_array.pkl'),'wb'))
    
    start = time.time()
    Results = match_wrapper(Fragments,Frag_Result,Seq_Result,out_path,ncpu=ncpu)
    Results.to_csv(os.path.join(out_path,'Results.csv'))
    print(round(time.time()-start,2),'Seconds to match fragments')
    
    #view_coverage(Results,Sequence,Colors = ['purple','blue','red'],step=30,shift=1)
    
    print('Total Run Time: ', round(time.time()-Master_Start,2))
    
    