#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 10:29:17 2019

Standardized plotting parameters for quality plots

@author: joshua
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib as mpl
from matplotlib.colors import ListedColormap

def style(ax, fontsize, labeltop=False, labelbottom=True, labelleft=True, labelright=False, c=None, linewidth=2):
    '''
    For plot/subplot frames
    '''
    
    # Check if this is a multipanel plot
    if hasattr(ax, '__len__'):
        ax = ax.ravel()
        for i in range(len(ax)):
            ax[i].xaxis.set_tick_params(which ='major', top=True, direction='in', width=2, labelsize=fontsize/1.5, length=8, labelbottom=labelbottom, labeltop=labeltop)
            ax[i].yaxis.set_tick_params(which ='major', right=True, direction='in', width=2, labelsize=fontsize/1.5, length=8, labelleft=labelleft, labelright=labelright)
            ax[i].xaxis.set_tick_params(which ='minor', top=True, direction='in', width=2, labelsize=fontsize/1.5, length=4, labelbottom=labelbottom, labeltop=labeltop)
            ax[i].yaxis.set_tick_params(which ='minor', right=True, direction='in', width=2, labelsize=fontsize/1.5, length=4, labelleft=labelleft, labelright=labelright)
            
            for axis in ['top','bottom','left','right']:
                ax[i].spines[axis].set_linewidth(2)
            
            if c is not None:
                for axis in ['top','bottom','left','right']:
                    ax[i].spines[axis].set_color(c)
                ax[i].xaxis.label.set_color(c)
                ax[i].yaxis.label.set_color(c)
                ax[i].tick_params(which='both', colors=c)
            
    else:
        ax.xaxis.set_tick_params(which ='major', top=True, direction='in', width=2, labelsize=fontsize/1.5, length=8, labelbottom=labelbottom, labeltop=labeltop)
        ax.yaxis.set_tick_params(which ='major', right=True, direction='in', width=2, labelsize=fontsize/1.5, length=8, labelleft=labelleft, labelright=labelright)
        ax.xaxis.set_tick_params(which ='minor', top=True, direction='in', width=2, labelsize=fontsize/1.5, length=4, labelbottom=labelbottom, labeltop=labeltop)
        ax.yaxis.set_tick_params(which ='minor', right=True, direction='in', width=2, labelsize=fontsize/1.5, length=4, labelleft=labelleft, labelright=labelright)

        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(linewidth)
        
        if c is not None:
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_color(c)
            ax.xaxis.label.set_color(c)
            ax.yaxis.label.set_color(c)
            ax.tick_params(which='both', colors=c)
            
        

def ticks(ax, xmajor, ymajor, xminor=None, yminor=None):
    '''
    Set major and minor tick spacing. The minor tick spacing is optional.
    The tick spacing may be entered as a single value or as an array of values
    for multiple subplots.
    '''
    if hasattr(ax, '__len__'):
        ax = ax.ravel()
        if not hasattr(xmajor, '__len__'):
            '''
            allow user to set single locator for all subplots.
            '''
            xmajor = [xmajor]*len(ax)
            ymajor = [ymajor]*len(ax)
        
        for i in range(len(ax)):
            ax[i].xaxis.set_major_locator(plt.MultipleLocator(xmajor[i]))
            ax[i].yaxis.set_major_locator(plt.MultipleLocator(ymajor[i]))
            if xminor is not None:
                if not hasattr(xminor, '__len__'):
                    xminor = [xminor]*len(ax)
                ax[i].xaxis.set_minor_locator(plt.MultipleLocator(xminor[i]))
            if yminor is not None:
                if not hasattr(yminor, '__len__'):
                    yminor = [yminor]*len(ax)
                ax[i].yaxis.set_minor_locator(plt.MultipleLocator(yminor[i]))
    else:       
        ax.xaxis.set_major_locator(plt.MultipleLocator(xmajor))
        ax.yaxis.set_major_locator(plt.MultipleLocator(ymajor))
        
        if xminor is not None:
            ax.xaxis.set_minor_locator(plt.MultipleLocator(xminor))
        if yminor is not None:
            ax.yaxis.set_minor_locator(plt.MultipleLocator(yminor))

def cbar_style(cbar, fontsize, side=False):
    '''
    For colorbars
    '''
    cbar.ax.tick_params(labelsize=fontsize/1.5, direction ='in', width=2)
    
    if cbar.orientation == 'horizontal':
        cbar.ax.tick_params(bottom=False, top=True, labelbottom=False, labeltop=True)
        if side is not False:
            cbar.ax.xaxis.set_label_position(side)
            cbar.ax.xaxis.set_ticks_position(side)
    elif cbar.orientation == 'vertical':
        cbar.ax.tick_params(left=False, right=True, labelleft=False, labelright=True)
        if side is not False:
            cbar.ax.yaxis.set_label_position(side)
            cbar.ax.yaxis.set_ticks_position(side)
    cbar.ax.set_frame_on(True)
    cbar.outline.set_linewidth(2)
    
    
def cbar_descrete(col_dict, labels):
    '''
    im - plt image to make cbar for
    fig - plt figure
    
    Returns
    cm - Custom Colormap
    fmt - colorbar format
    ticks - colorbar tick locations
    
    Example:
        
        col_dict={0:"navy",
              1:"cyan",
              2:"g",
              3:"orange",
              4:"r"}
        
        labels = ["Ret","SF","Comp","LINER", "Sey"]
        
        cm, fmt, ticks = cbar_descrete(col_dict, labels)
        fig.colorbar(im, format=fmt, ticks=ticks)
        
    '''    
    # We create a colorbar from our list of colors
    cm = ListedColormap([col_dict[x] for x in col_dict.keys()])
    # Let's also define the description of each category : 1 (blue) is Sea; 2 (red) is burnt, etc... Order should be respected here ! Or using another dict maybe could help.
    len_lab = len(labels)
    # Prepare bins for the normalizer
    norm_bins = np.sort([*col_dict.keys()]) + 0.5
    norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
    ## Make normalizer and formatter
    norm = mpl.colors.BoundaryNorm(norm_bins, len_lab, clip=True)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])
    # Give tick locations
    diff = norm_bins[1:] - norm_bins[:-1]
    tickz = norm_bins[:-1] + diff / 2
    
    #cb = fig.colorbar(im, format=fmt, ticks=tickz)
    return cm, fmt, tickz
    
def gaussian_den(x, y):
    '''
    Make a normalized number density cmap for scatter plots.
    '''
    mask = (~np.isnan(x))&(np.isfinite(x))&(~np.isnan(y))&(np.isfinite(y))
    
    x = x[mask]
    y = y[mask]
    
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    z = z/z.max() 
    return z

def legend(axs, fontsize=False):
    if fontsize==False:
        leg = axs.legend(fancybox=False, edgecolor='k')
    else:
        leg = axs.legend(fontsize=fontsize/2, fancybox=False, edgecolor='k')
    leg.get_frame().set_linewidth(2.0)