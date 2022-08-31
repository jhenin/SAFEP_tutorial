from alchemlyb.visualisation.dF_state import plot_dF_state
from alchemlyb.visualisation import plot_convergence
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plotGeneral(cumulative, perWindow, RT, width=8, height=4, PDFtype='KDE'):
    fig, ((cumAx, del1),( eachAx, del2),(hystAx, pdfAx)) = plt.subplots(3,2, sharex='col', sharey='row', gridspec_kw={'width_ratios': [2, 1]})

    fig.delaxes(del1)
    fig.delaxes(del2)

    # Cumulative change in kcal/mol
    cumAx.errorbar(cumulative.index, cumulative.BAR.f*RT, yerr=cumulative.BAR.errors, marker=None, linewidth=1)
    cumAx.set(ylabel=r'Cumulative $\rm\Delta G_{\lambda}$'+'\n(kcal/mol)')

    # Per-window change in kcal/mol
    eachAx.errorbar(perWindow.index, perWindow.BAR.df*RT, yerr=perWindow.BAR.ddf, marker=None, linewidth=1)
    eachAx.plot(perWindow.index, perWindow.EXP.dG_f*RT, marker=None, linewidth=1, alpha=0.5)
    eachAx.errorbar(perWindow.index, -perWindow.EXP.dG_b*RT, marker=None, linewidth=1, alpha=0.5)
    eachAx.set(ylabel=r'$\rm\Delta G_{\lambda}$'+'\n(kcal/mol)')

    
    #Hysteresis Plots
    diff = perWindow.EXP['difference']
    hystAx.vlines(perWindow.index, np.zeros(len(perWindow)), diff, label="fwd - bwd", linewidth=2)
    hystAx.set(xlabel=r'$\lambda$', ylabel=r'$\delta_\lambda$ (kcal/mol)', ylim=(-1,1))
    

    
    if PDFtype=='KDE':
        kernel = sp.stats.gaussian_kde(diff)
        pdfX = np.linspace(-1, 1, 1000)
        pdfY = kernel(pdfX)
        pdfAx.plot(pdfY, pdfX, label='KDE')
    elif PDFtype=='Histogram':
        pdfY, pdfX = np.histogram(diff, density=True)
        pdfX = pdfX[:-1]+(pdfX[1]-pdfX[0])/2
        pdfAx.plot(pdfY, pdfX,  label="Estimated Distribution")
    else:
        raise(f"Error: PDFtype {PDFtype} not recognized")
    
    pdfAx.set(xlabel=PDFtype)

    std = np.std(diff)
    mean = np.average(diff)
    temp = pd.Series(pdfY, index=pdfX)
    mode = temp.idxmax()
    
    textstr = r"$\rm mode=$"+f"{np.round(mode,2)}"+"\n"+fr"$\mu$={np.round(mean,2)}"+"\n"+fr"$\sigma$={np.round(std,2)}"
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    pdfAx.text(0.15, 0.95, textstr, transform=pdfAx.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    fig.set_figwidth(width)
    fig.set_figheight(height*3)
    fig.tight_layout()
    
    return fig, [cumAx,eachAx,hystAx,pdfAx] 



def doConvPlot(ax, X, fs, ferr, fwdColor, label=None):
    ax.errorbar(X, fs, yerr=ferr, marker=None, linewidth=1, color=fwdColor, markerfacecolor='white', markeredgewidth=1, markeredgecolor=fwdColor, ms=5, label=label)
    return ax



def convergencePlot(theax, fs, ferr, bs, berr, fwdColor='#0072B2', bwdColor='#D55E00', lgndF=None, lgndB=None):
    if not lgndF:
        lgndF=fwdColor
        lgndB=bwdColor
        
        
    lower = fs[-1]-ferr[-1]
    upper = fs[-1]+ferr[-1]
    theax.fill_between([0,1],[lower, lower], [upper, upper], color=bwdColor, alpha=0.25)
    theax.errorbar(np.arange(len(fs))/len(fs)+0.1, fs, yerr=ferr, marker='o', linewidth=1, color=fwdColor, markerfacecolor='white', markeredgewidth=1, markeredgecolor=fwdColor, ms=5)
    theax.errorbar(np.arange(len(bs))/len(fs)+0.1, bs, yerr=berr, marker='o', linewidth=1, color=bwdColor, markerfacecolor='white', markeredgewidth=1, markeredgecolor=bwdColor, ms=5, linestyle='--')


    
    theax.xaxis.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    
    finalMean = fs[-1]
    theax.axhline(y= finalMean, linestyle='-.', color='gray')
    theax.plot(0, finalMean, linewidth=1, color=lgndF, label='Forward Time Sampling')
    theax.plot(0, finalMean, linewidth=1, color=lgndB, linestyle='--', label='Backward Time Sampling')
    theax.set(xlabel='Fraction of Simulation Time', ylabel=r'Total $\rm\Delta G_{\lambda}$ (kcal/mol)')
    theax.legend()
    return theax



#Cannonical convergence plot
def convergence_plot(u_nk, tau=1, units='kT', RT=0.59):
    forward, forward_error, backward, backward_error = doConvergence(u_nk, num_points=10)

    if units=='kcal/mol':
        forward = forward*RT
        forward_error = forward_error*RT
        backward = backward*RT
        backward_error = backward_error*RT
    
    ax = plot_convergence(forward, forward_error, backward, backward_error)
    
    if units=='kcal/mol':
        ax.set(ylabel=r'$\rm\Delta G$'+'\n(kcal/mol)')

    return plt.gca()


def fb_discrepancy_plot(l_mid, dG_f, dG_b):
    plt.vlines(l_mid, np.zeros(len(l_mid)), dG_f + np.array(dG_b), label="fwd - bwd", linewidth=3)
    plt.legend()
    plt.title('Fwd-bwd discrepancies by lambda')
    plt.xlabel('Lambda')
    plt.ylabel('Diff. in delta-G')
    #Return the figure
    return plt.gca() 

def fb_discrepancy_hist(dG_f, dG_b):
    plt.hist(dG_f + np.array(dG_b));
    plt.title('Distribution of fwd-bwd discrepancies')
    plt.xlabel('Difference in delta-G')
    plt.ylabel('Count')
    #return the figure
    return plt.gca()
