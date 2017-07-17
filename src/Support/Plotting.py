"""!
@file Support/Plotting.py
@package Support

@defgroup Plotting Plotting

@brief A module containing general plotting functions

@author James Bevins

@date 14Jul17
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from DataAnalysis.Stats import red_chisq

#------------------------------------------------------------------------------#
def plot(*args, **kwargs):
    """!
    Plots a series of data with or without line and data points. The default
    behavior is to plot both points and markers. Data points with error bars
    are automatically plotted as points without lines.

    @param args: \e array \n
        Can provide up to 4 arrays of data in the format [x, y, uncert] to
        plot. \n
    @param kwargs: <em> optional plotting inputs </em> \n
        An optional list of additional plot options.  This is wrapped in
        kwargs because 2.7 doesn't support args and keyword specified
        arguments.  The options are listed as kwargs parameters below. \n
    @param includeLines: <em> kwargs boolean </em> \n
        Flag to include lines in the plot. By default, lines are included. \n
    @param includeMarkers: <em> kwargs boolean </em> \n
        Flag to use a markers on the plot. By default, markers are not
        included. \n
    @param logX: <em> kwargs boolean </em> \n
        Flag to use a log scale on the x axis \n
    @param logY: <em> kwargs boolean </em> \n
        Flag to use a log scale on the y axis \n
    @param title: <em> kwargs string </em> \n
        An optional specification for the plot title. \n
    @param xLabel: <em> kwargs string </em> \n
        An optional specification for the x axis label. \n
    @param yLabel: <em> kwargs string </em> \n
        An optional specification for the y axis label. \n
    @param dataLabel: <em> kwargs string </em> \n
        An optional specification for the data set label to use in the
        legend. \n
    @param legend: <em> kwargs string </em> \n
        An optional specification to include or not include a legend. \n
    @param legendLoc: <em> kwargs int </em> \n
        Specifies the legend location. \n
    @param savePath: <em> kwargs string </em> \n
        An optional specification for the save location. \n
    @param xMin: <em> kwargs integer or float </em> \n
        An optional specification for the minimum X axis value. \n
    @param xMax: <em> kwargs integer or float </em> \n
        An optional specification for the maximum X axis value. \n
    @param yMin: <em> kwargs integer or float </em> \n
        An optional specification for the minimum Y axis value. \n
    @param yMax: <em> kwargs integer or float </em> \n
        An optional specification for the maximum Y axis value. \n
    @param figsize: <em> kwargs tuple </em> \n
        The (x,y) scale of the plot. \n
    @param grid: <em> kwargs boolean </em> \n
        Specifies whether to plot the tick grid. \n
    @param xMinorTicks: <em> kwargs integer or float </em> \n
        An optional specification for the location of the X-axis minor tick
        marks. \n
    @param yMinorTicks: <em> kwargs integer or float </em> \n
        An optional specification for the location of the Y-axis minor tick
        marks. \n
    @param color: <em> kwargs list </em> \n
        Specifies the color cycle. \n
    @param linewidth: <em> kwargs list </em> \n
        Specifies the lineWidth cycle. \n
    @param linestyle: <em> kwargs list </em> \n
        Specifies the linestyle cycle. \n
    @param dashes: <em> kwargs list </em> \n
        Specifies the dashes cycle. \n
    @param marker: <em> kwargs list </em> \n
        Specifies the marker cycle. \n
    """

    # Set defaults if not specified since 2.7 sucks
    includeLines = kwargs.pop('includeLines', True)
    includeMarkers = kwargs.pop('includeMarkers', True)
    logX = kwargs.pop('logX', False)
    logY = kwargs.pop('logY', False)
    title = kwargs.pop('title', '')
    xLabel = kwargs.pop('xLabel', '')
    yLabel = kwargs.pop('yLabel', '')
    dataLabel = kwargs.pop('dataLabel', ['Data Set \#{}'.format(i) \
                                         for i in range(0, len(args))])
    legend = kwargs.pop('legend', True)
    legendLoc = kwargs.pop('legendLoc', 1)
    savePath = kwargs.pop('savePath', '')
    xMin = kwargs.pop('xMin', min(min(a[0]) for a in args)*0.75)
    xMax = kwargs.pop('xMax', max(max(a[0]) for a in args)*1.25)
    yMin = kwargs.pop('yMin', min(min(a[1]) for a in args)*0.5)
    yMax = kwargs.pop('yMax', max(max(a[1]) for a in args)*1.5)
    figsize = kwargs.pop('figsize', (9, 6))
    grid = kwargs.pop('grid', True)
    xMinorTicks = kwargs.pop('xMinorTicks', 0)
    yMinorTicks = kwargs.pop('yMinorTicks', 0)
    color = kwargs.pop('color', ['k', 'k', 'k', 'k', 'k', 'k'])
    linestyle = kwargs.pop('linestyle', ['-', ':', '-.', '--', '-', '-'])
    linewidth = kwargs.pop('linewidth', [2, 4])
    dashes = kwargs.pop('dashes', [[10, 0.1], [2, 2, 2, 2], [10, 5, 2, 5],
                                   [10, 5, 10, 5], [10, 2, 2, 2, 2, 2],
                                   [10, 2, 10, 2, 2, 2, 2, 2]])
    if includeMarkers:
        marker = kwargs.pop('marker', ['o', '^', '+', 's', 'd', '*', '>'])
    else:
        marker = [None]

    # Allow use of Tex sybols
    plt.rc('text', usetex=True)

    # Set up figure
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Set axes
    ax1.axis([xMin, xMax, yMin, yMax])
    if logX:
        ax1.set_xscale('log')
    if logY:
        ax1.set_yscale('log')

    # Set axes labels and plot title.
    ax1.set_title('{}'.format(title), fontsize=30, weight="bold")
    ax1.set_xlabel('{}'.format(xLabel), fontsize=22, weight='bold', y=-0.04)
    ax1.set_ylabel('{}'.format(yLabel), fontsize=22, weight="bold", x=-0.04)
    if xMinorTicks != 0:
        xMinorLocator = MultipleLocator(kwargs['xMinorTicks'])
        ax1.xaxis.set_minor_locator(xMinorLocator)
    if yMinorTicks != 0:
        yMinorLocator = MultipleLocator(kwargs['yMinorTicks'])
        ax1.yaxis.set_minor_locator(yMinorLocator)
    ax1.tick_params(axis='both', which='major', direction='inout',
                    labelsize=18, width=3, top=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', width=2.5,
                    top=True, right=True)
    ax1.tick_params(axis='both', which='major', direction='inout',
                    labelsize=18, width=3, top=True, right=True, length=7)
    ax1.tick_params(axis='both', which='minor', direction='in', width=2.5,
                   top=True, right=True, length=3)

    # Set up gridlines
    if grid:
        ax1.xaxis.grid(b=True, which='major', color='0.2', linestyle='-',
                       alpha=0.5)
        ax1.yaxis.grid(b=True, which='both', color='0.2', linestyle='-',
                       alpha=0.5)

    # Add datasets to plot
    n = 0
    l = 0
    for arg in args:
        if len(arg) == 3:
            if n >= len(dataLabel):
                ax1.errorbar(arg[0], arg[1], yerr=arg[2],
                             marker=marker[n%len(marker)], linestyle='None',
                             capsize=4, capthick=1.5, color=color[n%len(color)])
            else:
                ax1.errorbar(arg[0], arg[1], yerr=arg[2],
                             marker=marker[n%len(marker)], linestyle='None',
                             capsize=4, capthick=1.5, label=dataLabel[n],
                             color=color[n%len(color)])
        elif includeLines:
            ax1.plot(arg[0], arg[1], linewidth=linewidth[l//len(linestyle)],
                     color=color[n%len(color)],
                     linestyle=linestyle[l%len(linestyle)],
                     dashes=dashes[l%len(dashes)], marker=marker[n%len(marker)],
                     label=dataLabel[n])
            l += 1
        else:
            ax1.plot(arg[0], arg[1], linewidth=linewidth[n//len(linestyle)],
                     color=color[n%len(color)], linestyle='None',
                     marker=marker[n%len(marker)], label=dataLabel[n])
        n += 1

    # Add and locate legend
    if legend == True:
        plt.legend(borderaxespad=0.75, loc=legendLoc, fontsize=16,
                   handlelength=5, borderpad=0.5, labelspacing=0.75,
                   fancybox=True, framealpha=0.5, numpoints=1)

    plt.show()

    if savePath != '':
        fig.savefig(savePath, bbox_inches='tight')

#------------------------------------------------------------------------------#
def comp_plot(x, dataY, dataUncert, modelY, includeChi2=True,
              freeParams=1, **kwargs):
    """!
    Comparies a series of experimental data points (plots as points
    without lines) to a model (plots as a line).  The residuals are
    displayed (abs(residuals) if loxY=True), and the title can show
    the reduced Chi2 if specified.

    It is assumed that the x locations of the experimental data and
    model match.

    @param x: <em> list or array of integers or floats </em> \n
        The experimental data. \n
    @param dataY: <em> list or array of integers or floats </em> \n
        The experimental data. \n
    @param dataUncert: <em> list or array of integers or floats </em> \n
        Experimental data 1\f$\sigma\f$ standard deviation. \
    @param modelY: <em> list or array of integers or floats </em> \n
        The model data. \n.
    @param includeChi2: \e boolean \n
        Optional specifier to include the reduced chi2 calculation on
        the plot. \n
    @param freeParams: \e integer \n
        The number of free parameters in the model. \n
    @param kwargs: <em> optional plotting inputs </em> \n
        An optional list of additional plot options.The options are
        listed as kwargs parameters below \n
    @param logX: <em> kwargs boolean </em> \n
        Flag to use a log scale on the x axis \n
    @param logY: <em> kwargs boolean </em> \n
        Flag to use a log scale on the y axis \n
    @param title: <em> kwargs string </em> \n
        An optional specification for the plot title. \n
    @param dataLabel: <em> kwargs string </em> \n
        An optional specification for the data set label to use in the
        legend. \n
    @param legend: <em> kwargs string </em> \n
        An optional specification to include or not include a legend. \n
    @param xLabel: <em> kwargs string </em> \n
        An optional specification for the x axis label. \n
    @param yLabel: <em> kwargs string </em> \n
        An optional specification for the y axis label. \n
    @param savePath: <em> kwargs string </em> \n
        An optional specification for the save location. \n
    @param xMin: <em> kwargs integer or float </em> \n
        An optional specification for the minimum X axis value. \n
    @param xMax: <em> kwargs integer or float </em> \n
        An optional specification for the maximum X axis value. \n
    @param yMin: <em> kwargs integer or float </em> \n
        An optional specification for the minimum Y axis value. \n
    @param yMax: <em> kwargs integer or float </em> \n
        An optional specification for the maximum Y axis value. \n
    @param xMinorTicks: <em> kwargs integer or float </em> \n
        An optional specification for the location of the X-axis minor tick
        marks. \n
    @param yMinorTicks: <em> kwargs integer or float </em> \n
        An optional specification for the location of the Y-axis minor tick
        marks. \n
    """

    # Set defaults if not specified
    if 'logX' not in kwargs.keys():
        kwargs['logX'] = False
    if 'logY' not in kwargs.keys():
        kwargs['logY'] = False

    # Calculate the residuals
    if kwargs['logY']:
        residuals = abs(dataY - modelY)
    else:
        residuals = (dataY - modelY)

    # Set defaults if not specified
    if 'title' not in kwargs.keys():
        kwargs['title'] = ''
    if 'dataLabel' not in kwargs.keys():
        kwargs['dataLabel'] = ['experiment', 'model']
    if 'legend' not in kwargs.keys():
        kwargs['legend'] = True
    if 'xLabel' not in kwargs.keys():
        kwargs['xLabel'] = ''
    if 'yLabel' not in kwargs.keys():
        kwargs['yLabel'] = ''
    if 'savePath' not in kwargs.keys():
        kwargs['savePath'] = ''
    if 'xMin' not in kwargs.keys():
        kwargs['xMin'] = min(x)
    if 'xMax' not in kwargs.keys():
        kwargs['xMax'] = max(x)
    if 'yMin' not in kwargs.keys():
        kwargs['yMin'] = min(min(dataY), min(modelY), min(residuals))
        if kwargs['yMin'] < 0:
            kwargs['yMin'] = kwargs['yMin']*1.25
        else:
            kwargs['yMin'] = kwargs['yMin']*0.75
    if 'yMax' not in kwargs.keys():
        kwargs['yMax'] = max(max(dataY), max(modelY))*1.25
    if 'xMinorTicks' not in kwargs.keys():
        kwargs['xMinorTicks'] = 0
    if 'yMinorTicks' not in kwargs.keys():
        kwargs['yMinorTicks'] = 0

    # Allow use of Tex sybols
    plt.rc('text', usetex=True)

    # Set up figure
    fig = plt.figure(figsize=(9, 6))
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Set plot title.
    if includeChi2:
        redChiSq = red_chisq(dataY, modelY, dataUncert,
                             freeParams=freeParams)
        ax1.set_title('{} $\chi^2$={:.2f}'.format(kwargs['title'], redChiSq),
                     fontsize=30, weight="bold")
    else:
        ax1.set_title('{}'.format(kwargs['title']), fontsize=30, weight="bold")

    # Set primary axes
    ax1.axis([kwargs['xMin'], kwargs['xMax'], kwargs['yMin'], kwargs['yMax']])
    if kwargs['logX']:
        ax1.set_xscale('log')
    if kwargs['logY']:
        ax1.set_yscale('log')

    # Set axes labels
    ax1.set_xlabel('{}'.format(kwargs['xLabel']), fontsize=22,
                   weight='bold', y=-0.04)
    ax1.set_ylabel('{}'.format(kwargs['yLabel']), fontsize=22,
                   weight="bold", x=-0.04)
    if kwargs['xMinorTicks'] != 0:
        xMinorLocator = MultipleLocator(kwargs['xMinorTicks'])
        ax1.xaxis.set_minor_locator(xMinorLocator)
    if kwargs['yMinorTicks'] != 0:
        yMinorLocator = MultipleLocator(kwargs['yMinorTicks'])
        ax1.yaxis.set_minor_locator(yMinorLocator)
    ax1.tick_params(axis='both', which='major', direction='inout',
                    labelsize=18, width=3, top=True, right=True, length=7)
    ax1.tick_params(axis='both', which='minor', direction='in', width=2.5,
                   top=True, right=True, length=3)

    # Add datasets to plot
    ax1.errorbar(x, dataY, yerr=dataUncert, marker='o',
                 linestyle='None', capsize=4, capthick=1.5,
                 label=kwargs['dataLabel'][0], color='k')
    ax1.plot(x, modelY, linestyle='-', label=kwargs['dataLabel'][1],
             color='k')
    ax1.plot(x, residuals, linestyle='--', label='Residuals', color='k')

    # Add and locate legend
    if kwargs['legend'] == True:
        plt.legend(borderaxespad=0.75, loc=1, fontsize=16, handlelength=5,
                   borderpad=0.5, labelspacing=0.75, fancybox=True,
                   framealpha=0.5, numpoints=1)

    plt.show()

    if kwargs['savePath'] != '':
        fig.savefig(kwargs['savePath'], bbox_inches='tight')
