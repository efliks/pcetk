#!/usr/bin/python

import glob


pageheader = """# Gnuplot input file
#set   terminal pdf dashed size 16.6, 23.4 font "Verdana, 10" linewidth 2.0
set   terminal pdf dashed size 8.3, 11.7 font "Verdana, 10" linewidth 1.5
set   output "%s"

set   xrange [0:14]
set   yrange [0:1]
set   xtics  2.0 scale 0.8 offset 0,   0.3
set   mxtics
set   ytics  0.2 scale 0.8 offset 0.3, 0
set   mytics
set   xlabel "pH"           offset 0.0, 0.8  font "Verdana, 10"
set   ylabel "Probability"  offset 1.8, 0.0  font "Verdana, 10"

# B-deprotonated chromophore
set   style line  1  lc rgb "black"    lw 2     lt 1        # parent
set   style line  2  lc rgb "red"      lw 2     lt 1        # IFP 1.0
set   style line  3  lc rgb "blue"     lw 2     lt 1        # IFP 1.1
set   style line  4  lc rgb "green"    lw 2     lt 1        # IFP 1.2
set   style line  5  lc rgb "cyan"     lw 2     lt 1        # IFP 1.3
set   style line  6  lc rgb "magenta"  lw 2     lt 1        # IFP 1.4
set   style line  7  lc rgb "yellow"   lw 2     lt 1        # IFP 2.0

# Protonated chromophore
set   style line  8  lc rgb "black"    lw 4     lt 2        # parent
set   style line  9  lc rgb "red"      lw 4     lt 2        # IFP 1.0
set   style line 10  lc rgb "blue"     lw 4     lt 2        # IFP 1.1
set   style line 11  lc rgb "green"    lw 4     lt 2        # IFP 1.2
set   style line 12  lc rgb "cyan"     lw 4     lt 2        # IFP 1.3
set   style line 13  lc rgb "magenta"  lw 4     lt 2        # IFP 1.4
set   style line 14  lc rgb "yellow"   lw 4     lt 2        # IFP 2.0

set   multiplot layout 3, 2

unset key
"""

pagefooter = """
unset multiplot
"""

plotheader   = """
%s
set   size 0.33, 0.20
set   origin %f, %f
set   title "%s" font "Verdana Bold, 10" offset 0, -0.5
"""

plotfooter = """
%s
"""

plotcommand = "plot \\\n"


#===========================================================
def writepage (counter, data):
    content = [pageheader % ("page%02d.pdf" % counter)]
    content.extend (data)
    content.append (pagefooter)
    f = open ("page%02d.gnuplot" % counter, "w")
    f.writelines (content)
    f.close ()


#===========================================================
proteinlabels = ("curves_analytic", "curves_gmct", ) 
nproteins = len (proteinlabels)

files = []
for directory in proteinlabels:
    files.extend (glob.glob ("%s/*dat" % directory))


proteins = {}
for f in files:
    protein, filename = f.split ("/")
    segment, site, label = filename.split ("_")
    name   = site[:3]
    serial = int (site[3:])
    label  = label.split (".")[0]

    if not proteins.has_key (protein):
        sites = {}
        proteins[protein] = sites

    sites = proteins[protein]
    if sites.has_key (site):
        sname, sserial, instances = sites[site]
    else:
        instances = {}
        sites[site] = [name, serial, instances]
    instances[label] = f


instances = []
protein_items = proteins.iteritems ()
for protein_name, sites in protein_items:

    site_items = sites.iteritems ()
    for site_id, (site_name, site_serial, site_instances) in site_items:

        instance_items = site_instances.iteritems ()
        for instance_label, instance_file in instance_items:

            entry = [site_serial, site_name, instance_label]
            if entry not in instances:
                instances.append (entry)
instances.sort ()


#===========================================================
for counter, label in enumerate (proteinlabels):
    title = label
    plotcommand = "%s\"%%s\" w l ls %d t \"%s\"" % (plotcommand, counter + 1, title)
    if counter < (nproteins - 1):
        plotcommand = "%s, \\\n" % plotcommand
    else:
        plotcommand = "%s \n" % plotcommand


output = []
x, y   = 0, 0
page   = 1

for site_serial, site_name, instance_label in instances:
    ox = x * 0.33
    oy = 0.75 - y * 0.20
    title = "%s%d (%s)" % (site_name, site_serial, instance_label)
    sitekey = "%s%d" % (site_name, site_serial)

    filenames = []
    for protein_label in proteinlabels:
        sites = proteins[protein_label]

        if sites.has_key (sitekey):
            site = sites[sitekey]
            sname, sserial, sinstances = site
            filename = sinstances[instance_label]
        else:
            filename = "x"
        filenames.append (filename) 

    if x < 1 and y < 1:
        plotbefore = "set   key"
        plotafter  = "\nunset key\n"
    else:
        plotbefore = ""
        plotafter  = ""

    output.append (plotheader % (plotbefore, ox, oy, title))
    output.append (plotcommand % tuple (filenames))
    output.append (plotfooter % plotafter)

    x += 1
    if x > 2:
        x = 0
        y += 1
        if y > 3:
            writepage (page, output)
            y = 0
            page += 1
            output = []

if output:
  writepage (page, output)
