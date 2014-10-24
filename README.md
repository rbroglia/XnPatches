[![Ready in backlog](https://badge.waffle.io/szaghi/xnpatches.png?label=ready&title=Ready)](https://waffle.io/szaghi/xnpatches)
[![In Progress](https://badge.waffle.io/szaghi/xnpatches.png?label=in%20progress&title=In%20Progress)](https://waffle.io/szaghi/xnpatches)
[![Open bugs](https://badge.waffle.io/szaghi/xnpatches.png?label=bug&title=Open%20Bugs)](https://waffle.io/szaghi/xnpatches)

# XnPatches
### <a name="top">Xnavis PostProcessor
XnPatches loads Xnavis mesh and solution files and produces post-processed plot files.

## <a name="toc">Table of Contents

* [Team Members](#team-members)
* [What is XnPatches?](#what)
* [Todos](#todos)
* [Requirements](#requirements)
* [Copyrights](#copyrights)
* [Download XnPatches](#download)
* [Compiling Instructions](#compile)
* [Usage](#usage)
  + [Main Help](#usage-help)
  + [Post-processing only mesh files](#usage-only-mesh)
  + [Post-processing mesh and solutions files](#usage-mesh-sol)
  + [Utilities](#utilities)
* [Version History](#versions)

## <a name="team-members"></a>Team Members
* Stefano Zaghi    <stefano.zaghi@gmail.com>
* Riccardo Broglia

Go to [Top](#top) or [Toc](#toc)
## <a name="what"></a>What is XnPatches?
XnPatches loads Xnavis mesh and solution files and produces post-processed plot files. It is one of the many post-processors of Xnavis code.

Go to [Top](#top) or [Toc](#toc)
## <a name="todos"></a>Todos
+ Refactor utilities;
+ any feature request is welcome!

Go to [Top](#top) or [Toc](#toc)
## <a name="download"></a>Download XnPatches
If you use `git` it could be convenient to clone this repository:
```bash
git clone https://github.com/szaghi/XnPatches
```
Other 2 possibilities are:

1. use the GitHub **Download ZIP** button on the right sidebar of this page;
2. download one of the releases in the [release page](https://github.com/szaghi/XnPatches/releases), also listed at the end of this page.

Go to [Top](#top) or [Toc](#toc) 
## <a name="requirements"></a>Requirements
+ Modern Fortran Compiler (standard 2003+);
+ a lot of patience with the author.

XnPatches is developed on a GNU/Linux architecture. For Windows architecture there is no support, however it should be work out-of-the-box.

Go to [Top](#top) or [Toc](#toc)
## <a name="Copyrights"></a>Copyrights
XnPatches is an open source project, it is distributed under the [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html). Anyone is interest to use, to develop or to contribute to XnPatches is welcome.

Go to [Top](#top) or [Toc](#toc)
## <a name="compile"></a>Compiling Instructions
XnPatches has been developed on GNU/Linux architectures. Other OS are not supported (and in general there is no best alternative to GNU/Linux :-).

XnPatches have been successfully compiled with the following compilers:

+ GNU gfortran (version 4.7.0 or higher);
+ Intel Fortran Compiler ifort (version 12.0 or higher)


XnPatches is constituted by several modules. Therefore there are many dependences. The most easy way to compile the code is to start with the provided makefile thus it is necessary that the system has "Make" program (preferably GNU make http://www.gnu.org/software/make/).

The provided makefile has several options. It has one rule that prints all options available and the default settings. Typing in the shell prompt: `code make help` the following output will be printed:

```bash
 Make options of XnPatches code

 Compiler choice: COMPILER=intel => default
  COMPILER=gnu   => GNU gfortran           
  COMPILER=intel => Intel Fortran         

 Compiling options
  DEBUG=yes(no)    => on(off) debug                  (default no)
  F03STD=yes(no)   => on(off) check standard fortran (default no)
  OPTIMIZE=yes(no) => on(off) optimization           (default no)
  OPENMP=yes(no)   => on(off) OpenMP directives      (default no)                                                                                                                                                                                                     
  BIGEIN=yes(no)   => on(off) Big Endian input files (default yes)                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                      
 Preprocessing options                                                                                                                                                                                                                                                
  R16P=yes(no) => on(off) definition of real with "128 bit" (default no)                                                                                                                                                                                              
                                                                                                                                                                                                                                                                      
 Executable directory                                                                                                                                                                                                                                                 
  DEXE="your_path" => directory where exe is placed (default ~/bin/)                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                      
 External libraries                                                                                                                                                                                                                                                   
  TECIO=yes(no) => on(off) Tecplot IO library linking (default )                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                      
 Provided Rules                                                                                                                                                                                                                                                       
  Defualt rule    => ~/bin/XnPatches                                                                                                                                                                                                                                  
  help            => printing this help message                                                                                                                                                                                                                       
  ~/bin/XnPatches => building OFF code                                                                                                                                                                                                                                
  cleanobj        => cleaning compiled object                                                                                                                                                                                                                         
  cleanmod        => cleaning .mod files                                                                                                                                                                                                                              
  cleanmsg        => cleaning make-log massage files                                                                                                                                                                                                                  
  cleanexe        => cleaning executable files                                                                                                                                                                                                                        
  clean           => running cleanobj, cleanmod and cleanmsg                                                                                                                                                                                                          
  cleanall        => running clean and cleanexe                                                                                                                                                                                                                       
  tar             => creating a tar archive of the project
```
For example compiling in debug mode with the Intel Fortran compiler you can type:
```bash
make DEBUG=yes COMPILER=intel
```

Go to [Top](#top) or [Toc](#toc)
## <a name="usage"></a>Usage
### <a name="usage-help">Main Help
XnPatches is is a Command Line Tool. To list the available options run it as following:
```bash
./XnPatches -h
```
this help will echoed
```bash
 XnPatches
 Post processing code for Xnavis code
 Usage:
   XnPatches -g file_geo -i file_icc
            [-o file_output
             -s file_solution
             -p #bc_of_processed_patch
             -cell
             -ls
             -eq #turbulent_eq_model(0,1,2)
             -Re #Re
             -Fr #Fr
             -zfs #zfs
             -forces
             -forcesRB
             -yplus
             -metrics
             -ascii
             -tec yes/no
             -vtk yes/no
             -proc #proc
             -os UIX/WIN]
 
 Optional arguments meaning and default values:
  -o file_output   => output file name; default is basename of grd file with the proper extension
  -s file_solution => solution file name; if passed the solution variables are saved
  -p 1             => walls are default processed patches
  -cell            => all variables other than nodes coord. are cell centered (default no, node centered)
  -ls              => solution with level set, Froude number and free surface quote are necessary if
                      forces are computed (default no)
  -eq 1
  -Re #Re          => Reynolds number
  -Fr #Fr          => Froude number
  -zfs #zfs        => Quote of free surface
  -forces          => compute forces, Reynolds number is necessary (default no, do not compute them)
  -forcesB         => compute forces and save also RB outputs (default no, do not compute them)
  -yplus           => compute y+, Reynolds number ans solution are necessary (default no, do not compute)
  -metrics         => save metrics variables (default no, do not save metrics)
  -ascii           => write ascii output file (default no, write binary one)
  -tec yes/no      => write (or not) Tecplot file format (default yes)
  -vtk yes/no      => write (or not) VTK file format (default no)
  -proc #proc      => if file "proc.input" if found global blocks numeration is used; #proc is the process
                      number of the current processed file
  -os UIX/WIN      => type of Operating System write (default *UIX OS type)
 
 Examples:
  XnPatches -g cc.01.grd -i cc.01                              -o mesh.01 (process only mesh of wall)
  XnPatches -g cc.01.grd -i cc.01 -s sol.00.01 -forces -Re 1d6 -o sol.01  (solution and forces are saved)
 
 Note:
   1) the output file name extension is not necessary because it assigned according to the type of output:
      binary       Tecplot => .plt
      ascii        Tecplot => .dat
      binary/ascii VTK     => .vtm
   2) if a file name "mb.par" is present into the directory where XnPatches is executed the values of Re,
      Fr, zfs and the turbulence model are loaded from this file thus they can be omitted from the command
      line arguments list.
   3) if forces are computed output file with suffix "-forces.dat" is also saved: this file contains the
      integral of the forces over all the patches and the "wet" surface.
   4) all the variables other than the nodes coordinates are saved at cell center if "-cell" option is used;
      if blanking is used the blanking mode must be "any corners" or "primary values".
```
### <a name="usage-only-mesh">Post-processing only mesh files
#### GRD files (no ghost cells)
Not yet supported.
#### Overset/Overott files (with ghost cells)
```bash
  XnPatches -g cc.01.grd -i cc.01 -o mesh
```
By default the wall boundary conditions, `1`, are extracted thus a file named `mesh.01.plt` is generated, where the suffix `.01` indicates the boundary conditions. 
```bash
  XnPatches -g cc.01.grd -i cc.01 -o inflow -p 3
```
the inflow `3` patches are extracted into the file `inflow.03.plt`
### <a name="usage-mesh-sol">Post-processing mesh and solutions files
```bash
  XnPatches -g cc.01.grd -i cc.01 -s sol.00.01 -o sol
```
A file named `sol.01.plt` is generated where the wall patches are extracted with also the solutions fields on that patches. 
### <a name="utilities">Utilities
To be written.

Go to [Top](#top) or [Toc](#toc)

## <a name="versions"></a>Version History
In the following the changelog of most important releases is reported.
### v0.0.2 
##### Download [ZIP](https://github.com/szaghi/XnPatches/archive/v0.0.2.zip) ball or [TAR](https://github.com/szaghi/XnPatches/archive/v0.0.2.tar.gz) one
XnPatches.f90 module split for improve compiling speed. Fully backward compatible.
### v0.0.1 
##### Download [ZIP](https://github.com/szaghi/XnPatches/archive/v0.0.1.zip) ball or [TAR](https://github.com/szaghi/XnPatches/archive/v0.0.1.tar.gz) one
Stable Release. Fully backward compatible.
