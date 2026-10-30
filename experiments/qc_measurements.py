"""Descriptive measurements only. No tissue diagnosis or acceptance thresholds."""
import numpy as np

def measure(arr):
 x=np.asarray(arr,dtype=np.float32);assert x.ndim==3 and x.shape[2]==3
 gray=.299*x[...,0]+.587*x[...,1]+.114*x[...,2]
 lap=-4*gray.copy();lap[1:]+=gray[:-1];lap[:-1]+=gray[1:];lap[:,1:]+=gray[:,:-1];lap[:,:-1]+=gray[:,1:]
 interior=lap[1:-1,1:-1]
 mx=x.max(2);mn=x.min(2);sat=np.divide(mx-mn,mx,out=np.zeros_like(mx),where=mx>0);value=mx/255
 old_tissue=((sat>.04)&(value<.985))|(value<.80)
 chromatic=(sat>.04)&(value<.985)
 # Fraction admitted solely by the darkness branch of the old tissue formula.
 dark_only=old_tissue&~chromatic
 return {'focus_zero_padded_var':float(lap.var()),'focus_interior_var':float(interior.var()),'gray_std':float(gray.std()),'mean_saturation':float(sat.mean()),'mean_channel_range':float((mx-mn).mean()),'old_tissue_fraction':float(old_tissue.mean()),'chromatic_fraction_same_old_cutoff':float(chromatic.mean()),'dark_only_fraction':float(dark_only.mean()),'interior_gradient_mean':float((abs(np.diff(gray[1:-1,1:-1],axis=0)).mean()+abs(np.diff(gray[1:-1,1:-1],axis=1)).mean())/2)}
