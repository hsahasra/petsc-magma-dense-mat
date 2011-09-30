#!/usr/bin/env python
import PETSc.package

class Configure(PETSc.package.NewPackage):
  def __init__(self, framework):
    PETSc.package.NewPackage.__init__(self, framework)
    self.download         = ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/ParMetis-3.2.0-p1.tar.gz']
    self.functions        = ['ParMETIS_V3_PartKway']
    self.includes         = ['parmetis.h']
    self.liblist          = [['libparmetis.a','libmetis.a']]
    self.needsMath        = 1
    self.double           = 0
    self.complex          = 1
    self.requires32bitint = 0 # see note in configureLibrary() method below
    self.worksonWindows    = 1
    self.downloadonWindows = 1
    return

  def setupDependencies(self, framework):
    PETSc.package.NewPackage.setupDependencies(self, framework)
    self.deps = [self.mpi]
    return

  def Install(self):
    import os
    import sys

    makeinc        = os.path.join(self.packageDir,'make.inc')
    installmakeinc = os.path.join(self.confDir,'ParMetis')
    metisconfigheader    = os.path.join(self.packageDir,'METISLib','configureheader.h')
    parmetisconfigheader = os.path.join(self.packageDir,'ParMETISLib','configureheader.h')
    
    # Configure ParMetis 
    g = open(makeinc,'w')
    g.write('SHELL          = '+self.programs.SHELL+'\n')
    g.write('CP             = '+self.programs.cp+'\n')
    g.write('RM             = '+self.programs.RM+'\n')
    g.write('MKDIR          = '+self.programs.mkdir+'\n')

    g.write('AR             = '+self.setCompilers.AR+'\n')
    g.write('ARFLAGS        = '+self.setCompilers.AR_FLAGS+'\n')
    g.write('AR_LIB_SUFFIX  = '+self.setCompilers.AR_LIB_SUFFIX+'\n')
    g.write('RANLIB         = '+self.setCompilers.RANLIB+'\n')

    g.write('PARMETIS_ROOT  = '+self.packageDir+'\n')
    g.write('PREFIX         = '+self.installDir+'\n')
    g.write('METISLIB       = $(PARMETIS_ROOT)/libmetis.$(AR_LIB_SUFFIX)\n')
    g.write('PARMETISLIB    = $(PARMETIS_ROOT)/libparmetis.$(AR_LIB_SUFFIX)\n')
    
    self.setCompilers.pushLanguage('C')
    cflags = self.setCompilers.getCompilerFlags().replace('-Wall','').replace('-Wshadow','')
    cflags += ' ' + self.headers.toString(self.mpi.include)+' '+self.headers.toString('.')
        
    g.write('CC             = '+self.setCompilers.getCompiler()+'\n')
    g.write('CFLAGS         = '+cflags+'\n')
    # parmetis uses defaut 'make' targets, and this uses TARGET_ARCH var. If this var
    # is set incorrectly in user env - build breaks.
    g.write('TARGET_ARCH    = \n')

    self.setCompilers.popLanguage()
    g.close()

#   Warning: PTScotch also installs a file metis.h which will be incorrectly found by ParMetis when compiling so you cannot build PETSc to use
#   both PTScotch and ParMetis in the same PETSc build
#    if self.framework.argDB['download-scotch']:
#      raise RuntimeError('Cannot use both --download-scotch and --download-parmetis')
   
    if self.installNeeded('make.inc'):    # Now compile & install
      self.framework.outputHeader(metisconfigheader,prefix='METIS')
      self.framework.outputHeader(parmetisconfigheader,prefix='PARMETIS')
      try:
        self.logPrintBox('Compiling & installing Parmetis; this may take several minutes')
        output,err,ret  = PETSc.package.NewPackage.executeShellCommand('cd '+self.packageDir+' && make clean && make lib && make minstall && make clean', timeout=2500, log = self.framework.log)
      except RuntimeError, e:
        raise RuntimeError('Error running make on ParMetis: '+str(e))
      self.postInstall(output+err,'make.inc')
    return self.installDir

  def configureLibrary(self):
    PETSc.package.NewPackage.configureLibrary(self)
    #  PETSc can use ParMetis directly if built with 32 bit indices; ParMetis cannot work with 64 bit indices
    #  so we turn off useddirectly, but it can still be used by SuperLU_Dist
    if self.libraryOptions.integerSize == 64:
      self.useddirectly = 0
