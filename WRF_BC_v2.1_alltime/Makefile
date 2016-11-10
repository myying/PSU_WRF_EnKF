#	Top-level Makefile for UPDATE_WRF_BC

#	Macros, these should be generic for all machines

.IGNORE:

AR	=	ar ru
CD	=	cd
LN	=	ln -s
MAKE	=	make -f Makefile
RM	=	/bin/rm -f
RM_LIST	=	*.o core *.i *.mod fort.* *.out namelist.* *~ *.exe

#	Targets for supported architectures

default:	make_rules
		find . -name \*.f -exec rm {} \;
		( $(CD) src ; $(MAKE) update_wrf_bc.exe )
		( $(RM) update_wrf_bc.exe ;   $(LN) src/update_wrf_bc.exe . )

clean:
	$(RM) $(RM_LIST)
	($(CD) src ; $(MAKE) clean )

make_rules:
	uname -a > .tmpfile
	grep OSF .tmpfile								; \
	if [ $$? = 0 ] ; then echo "Compiling for Compaq"				; \
		echo "CPP	= /usr/bin/cpp"			>  compiler_macros	; \
		echo "CPPFLAGS	= -I. -C -P -DDEC"		>> compiler_macros	; \
		echo "FC	= f90"				>> compiler_macros	; \
		echo "FCFLAGS	= -C -free -O2 -fpe0 -convert big_endian -pg -g1" >> compiler_macros	; \
		echo "LDFLAGS	= -O4 -pg"			>> compiler_macros	; \
		echo "FPPFLAG	= -DDEC -I/usr/local/netcdf/include -I./"			>> compiler_macros	; \
		echo "LIBS	= -L/usr/local/netcdf/lib -lnetcdf -lnetcdff"	>> compiler_macros	; \
	else grep AIX .tmpfile	 							; \
	if [ $$? = 0 ] ; then echo "Compiling for IBM"						; \
		echo "CPP	= /usr/lib/cpp"			>  compiler_macros	; \
		echo "CPPFLAGS	= -I. -C -P -DIBM"		>> compiler_macros	; \
		echo "FC	= xlf90"			>> compiler_macros	; \
		echo "FCFLAGS	= -qlanglvl=90pure -O3 -qarch=auto -qnosave -qmaxmem=-1 -Q"	>> compiler_macros	; \
		echo "LDFLAGS	= -O2" 				>> compiler_macros	; \
		echo "FPPFLAG	= -DIBM -I/usr/local/include -I."	>> compiler_macros	; \
		echo "LIBS	= -L/usr/local/lib32/r4i4 -lnetcdf -lnetcdff"	>> compiler_macros	; \
	else grep Darwin .tmpfile	 							; \
	if [ $$? = 0 ] ; then echo "Compiling for Mac OS X"						; \
		echo "CPP	= /usr/bin/cpp"			>  compiler_macros	; \
		echo "CPPFLAGS	= -I. -C -P -DMAC"		>> compiler_macros	; \
		echo "FC	= xlf90"			>> compiler_macros	; \
		echo "FCFLAGS	= -qlanglvl=90pure -O3 -qarch=auto -qnosave -qmaxmem=-1 -Q"	>> compiler_macros	; \
		echo "LDFLAGS	= -O2" 				>> compiler_macros	; \
		echo "FPPFLAG	= -DMAC -I/usr/local/netcdf-3.5.1-xlf/include -I."	>> compiler_macros	; \
		echo "LIBS	= -L/usr/local/netcdf-3.5.1-xlf/lib -lnetcdf -lnetcdff"	>> compiler_macros	; \
	else grep Linux .tmpfile 							; \
	if [ $$? = 0 ] ; then echo "Compiling for Linux"				; \
		echo "CPP	= /lib/cpp"			>  compiler_macros	; \
		echo "CPPFLAGS	= -I. -C -P -DLINUX -traditional -Dlinux"	>> compiler_macros	; \
		echo "FC	= ifort"			>> compiler_macros	; \
		echo "FCFLAGS	= -FR -convert big_endian"	>> compiler_macros	; \
		echo "LDFLAGS	= " 				>> compiler_macros	; \
		echo "FPPFLAG	= -DLINUX -I$(NETCDF)/include -I."	>> compiler_macros	; \
		echo "LIBS	= -L$(NETCDF)/lib -lnetcdf -lnetcdff"	>> compiler_macros	; \
	else echo "Do not know how to compile for the `cat .tmpfile` machine." 		; \
	fi ; \
	fi ; \
	fi ; \
	fi

	echo "AR	= $(AR)"	>> compiler_macros
	echo "RM	= $(RM)"	>> compiler_macros
	echo "RM_LIST	= $(RM_LIST)"	>> compiler_macros
	echo "CD	= $(CD)"	>> compiler_macros
	echo "LN	= $(LN)"	>> compiler_macros
	echo "MAKE	= $(MAKE)"	>> compiler_macros
	echo "SHELL	= /bin/sh"	>> compiler_macros
	echo "TOUCH	= touch"	>> compiler_macros

