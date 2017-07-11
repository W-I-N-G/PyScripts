"""!
@file Support/Plotting.py
@package Support

@defgroup Plotting Plotting

@brief A module containing general plotting functions

@author James Bevins

@date 11Jul17
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#------------------------------------------------------------------------------#
def line_plot(*args, **kwargs):
    """!
    Plots a series of lines without data points.

    @param args: \e array \n
        Can provide up to 4 arrays of data in the format [x, y, uncert] to
        plot. \n
    @param kwargs: <em> optional plotting inputs </em> \n
        An optional list of additional plot options.  This is wrapped in
        kwargs because 2.7 doesn't support args and keyword specified
        arguments.  The options are listed as kwargs parameters below \n
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

    # Set defaults if not specified since 2.7 sucks
    if 'logX' not in kwargs.keys():
        kwargs['logX'] = False
    if 'logY' not in kwargs.keys():
        kwargs['logY'] = False
    if 'title' not in kwargs.keys():
        kwargs['title'] = ''
    if 'dataLabel' not in kwargs.keys():
        kwargs['dataLabel'] = []
        for i in range(0, len(args)):
            kwargs['dataLabel'].append('Data Set \#{}'.format(i))
    if 'legend' not in kwargs.keys():
        kwargs['legend'] = True
    if 'xLabel' not in kwargs.keys():
        kwargs['xLabel'] = ''
    if 'yLabel' not in kwargs.keys():
        kwargs['yLabel'] = ''
    if 'savePath' not in kwargs.keys():
        kwargs['savePath'] = ''
    if 'xMin' not in kwargs.keys():
        kwargs['xMin'] = min(min(a[0] for a in args))*0.75
    if 'xMax' not in kwargs.keys():
        kwargs['xMax'] = max(max(a[0] for a in args))*1.25
    if 'yMin' not in kwargs.keys():
        kwargs['yMin'] = min(min(a[1] for a in args))*0.5
    if 'yMax' not in kwargs.keys():
        kwargs['yMax'] = max(max(a[1] for a in args))*1.5
    if 'xMinorTicks' not in kwargs.keys():
        kwargs['xMinorTicks'] = 0
    if 'yMinorTicks' not in kwargs.keys():
        kwargs['yMinorTicks'] = 0

    # Allow use of Tex sybols
    plt.rc('text', usetex=True)

    # Set up figure
    fig = plt.figure(figsize=(9, 6))
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Preset data set format scheme
    linewidth = [2.5]
    linestyle = ['-', '-', ':', '-.']
    dashes = [[10, 0.1], [5, 2, 5, 2], [0.5, 0.5], [10, 2.5, 1, 2.5]]
    color = ['k', 'k', 'k', 'k']

    # Set axes
    ax1.axis([kwargs['xMin'], kwargs['xMax'], kwargs['yMin'],
              kwargs['yMax']])
    if kwargs['logX']:
        ax1.set_xscale('log')
    if kwargs['logY']:
        ax1.set_yscale('log')

    # Set axes labels and plot title.
    ax1.set_title('{}'.format(kwargs['title']), fontsize=30, weight="bold")
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
                    labelsize=18, width=3, top=True, right=True)
    ax1.tick_params(axis='both', which='minor', direction='in', width=2.5,
                    top=True, right=True)
    ax1.tick_params(axis='both', which='major', direction='inout',
                    labelsize=18, width=3, top=True, right=True, length=7)
    ax1.tick_params(axis='both', which='minor', direction='in', width=2.5,
                   top=True, right=True, length=3)

    # Add datasets to plot
    n = 0
    for arg in args:
        ax1.plot(arg[0], arg[1], linewidth=linewidth[0], color=color[n],
                 linestyle=linestyle[n], dashes=dashes[n], marker=None,
                 label=kwargs['dataLabel'][n])
        if len(arg) == 3:
            ax1.errorbar(arg[0], arg[1], yerr=arg[2], marker=None,
                         linestyle='None', capsize=4, capthick=1.5)
        n += 1

    # Add and locate legend
    if kwargs['legend'] == True:
        plt.legend(borderaxespad=0.75, loc=1, fontsize=16, handlelength=5,
                   borderpad=0.5, labelspacing=0.75, fancybox=True,
                   framealpha=0.5, numpoints=1)

    plt.show()

    if kwargs['savePath'] != '':
        fig.savefig(kwargs['savePath'], bbox_inches='tight')

#------------------------------------------------------------------------------#
def scatter_plot(*args, **kwargs):
    """!
    Plots a series of data points without lines.

    @param args: \e array \n
        Can provide up to 7 arrays of data in the format [x, y, uncert] to
        plot. \n
    @param kwargs: <em> optional plotting inputs </em> \n
        An optional list of additional plot options.  This is wrapped in
        kwargs because 2.7 doesn't support args and keyword specified
        arguments.  The options are listed as kwargs parameters below \n
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

    # Set defaults if not specified since 2.7 sucks
    if 'logX' not in kwargs.keys():
        kwargs['logX'] = False
    if 'logY' not in kwargs.keys():
        kwargs['logY'] = False
    if 'title' not in kwargs.keys():
        kwargs['title'] = ''
    if 'dataLabel' not in kwargs.keys():
        kwargs['dataLabel'] = []
        for i in range(0, len(args)):
            kwargs['dataLabel'].append('Data Set \#{}'.format(i))
    if 'legend' not in kwargs.keys():
        kwargs['legend'] = True
    if 'xLabel' not in kwargs.keys():
        kwargs['xLabel'] = ''
    if 'yLabel' not in kwargs.keys():
        kwargs['yLabel'] = ''
    if 'savePath' not in kwargs.keys():
        kwargs['savePath'] = ''
    if 'xMin' not in kwargs.keys():
        kwargs['xMin'] = min(min(a[0] for a in args))*0.75
    if 'xMax' not in kwargs.keys():
        kwargs['xMax'] = max(max(a[0] for a in args))*1.25
    if 'yMin' not in kwargs.keys():
        kwargs['yMin'] = min(min(a[1] for a in args))*0.5
    if 'yMax' not in kwargs.keys():
        kwargs['yMax'] = max(max(a[1] for a in args))*1.5
    if 'xMinorTicks' not in kwargs.keys():
        kwargs['xMinorTicks'] = 0
    if 'yMinorTicks' not in kwargs.keys():
        kwargs['yMinorTicks'] = 0

    # Allow use of Tex sybols
    plt.rc('text', usetex=True)

    # Set up figure
    fig = plt.figure(figsize=(9, 6))
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Preset data set format scheme
    marker = ['o', '^', '+', 's', 'd', '*', '>']

    # Set axes
    ax1.axis([kwargs['xMin'], kwargs['xMax'], kwargs['yMin'],
              kwargs['yMax']])
    if kwargs['logX']:
        ax1.set_xscale('log')
    if kwargs['logY']:
        ax1.set_yscale('log')

    # Set axes labels and plot title.
    ax1.set_title('{}'.format(kwargs['title']), fontsize=30, weight="bold")
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
    n = 0
    for arg in args:
        if len(arg) == 3:
            ax1.errorbar(arg[0], arg[1], yerr=arg[2], marker=marker[n],
                         linestyle='None', capsize=4, capthick=1.5,
                         label=kwargs['dataLabel'][n], color='k')
        else:
            ax1.plot(arg[0], arg[1], marker[n], linestyle='None',
                     label=kwargs['dataLabel'][n], color='k')
        n += 1

    # Add and locate legend
    if kwargs['legend'] == True:
        plt.legend(borderaxespad=0.75, loc=1, fontsize=16, handlelength=5,
                   borderpad=0.5, labelspacing=0.75, fancybox=True,
                   framealpha=0.5, numpoints=1)

    plt.show()

    if kwargs['savePath'] != '':
        fig.savefig(kwargs['savePath'], bbox_inches='tight')
