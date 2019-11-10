#!/usr/bin/env python2

import sys,os,platform
import publishfunctions

curdir = os.path.dirname(os.path.realpath(__file__))

# ------- Configuration parameters -------------

projectname='gspbox'

if platform.system()=='Darwin':
    homefolder='/Users/user'
else:
    homefolder='/home/user'

project=homefolder+'/'+projectname+'/'


# -------- Configuration of mat2doc ------------
mat2docpath='~/mat2doc'

try :
    from publish_local import *
except:
    pass

# -------- Automatique configuration ------------
import conf
outputdir=conf.outputdir
outputdirphp=outputdir+projectname+'-php/'
outputdirmat=outputdir+projectname+'-mat/'
outputdirrelease=outputdir+projectname+'-release/'
outputdirtex=outputdir+projectname+'-tex/'
outputdirhtml=outputdir+projectname+'-html/'



f=file(project+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

# ------- do not edit below this line ----------
f = open(project + 'mat2doc/startup.m', 'w')
f.write('addpath ' + unlocxpath + '\n')
f.write('init_unlocbox;\n\n')
f.write('addpath ' + gsppath + '\n')
f.write('gsp_start;\n');
f.close()


todo=sys.argv

if 'package' in todo:
    todo.append('release')

#  Optional parameters

if 'fast' in todo:
    plot='--no-execplot '    
else:
    plot='--execplot '


if 'rebuild' in todo:
    build='--rebuild '
elif 'cached' in todo:
    build='--cached '
else:
    build='--auto '



#  Publish
for mode in  ['mat', 'php', 'html', 'tex']:
    if mode in todo:
        s = '%s %s/mat2doc.py %s%s %s %s' % ('PYTHONPATH="%s:$PYTHONPATH"' % (curdir,), mat2docpath, plot if mode != 'mat' else '', build if mode != 'mat' else '', project, mode,)

        os.system(s)

if 'tex' in todo:
    s = 'cp '+ project + 'mat2doc/tex/main/* '+outputdirtex
    os.system(s)
    s = 'cp '+ project + 'mat2doc/project.bib '+outputdirtex
    os.system(s)
    s = 'cd '+ outputdirtex +' && make clean && make'
    os.system(s)


# Release version
if 'release' in todo or 'package' in todo or 'pushrelease' in todo:
    os.system('rm -r '+outputdirrelease)
    os.system('cp -r '+outputdirmat+' '+outputdirrelease)

    os.system('rm -r '+outputdirrelease+'private')
    os.system('rm -r '+outputdirrelease+'oose')
    os.system('rm -r '+outputdirrelease+'testing')
    # os.system('rm -r '+outputdirrelease+'graph_ml')
    os.system('rm -r '+outputdirrelease+'test_gsptoolbox')
    os.system('rm -r '+outputdirrelease+'to\ be\ removed')
    os.system('rm -r '+outputdirrelease+'.git')
    os.system('rm -r '+outputdirrelease+'.gitignore')


#  Packaging
if 'package' in todo:
    
    fname=outputdir+projectname+'-'+versionstring
    # Create the Unix src package
    os.system('rm '+fname+'.tar.gz')
    os.system('cd '+outputdir+' && cp -r ' +projectname+'-release '+projectname+' && tar zcvf '+fname+'.tar.gz '+projectname+'/')
   
    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    publishfunctions.unix2dos(outputdir + projectname +'/')
    os.system('cd '+outputdir+' && zip -r '+fname+'.zip '+projectname+'/')    

    os.system('rm -r '+outputdir+projectname)

if 'copyhtml' in todo:
    htmldirgit = homefolder+'/work/git/website/gspbox-html/'
    # os.system('cp -r '+outputdirhtml+'include/* '+ htmldirgit+'include/')
    # os.system('rm -r '+outputdirhtml+'include ')
    os.system('cp -r '+outputdirhtml+'* '+ htmldirgit+'doc/')
    print('cp -r '+outputdirhtml+'* '+ htmldirgit+'doc/')

#  Send to the server

# if 'pushrelease' in todo:
#     s='rsync -r --progress '+outputdirrelease+'* '+outputdircode
#     print s    
#     os.system(s)  

# if 'pushdoc' in todo:
#     s='rsync -av '+outputdirphp+'* '+outputdirweb+'doc/'
#     print s    
#     os.system(s)  








