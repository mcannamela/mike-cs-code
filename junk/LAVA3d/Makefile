SRCS = \
	BCCC.f BCMOM.f CARRY.f CCFLUX.f CCFLUY.f CCFLUZ.f CHECK.f CHEM.f \
	CHEMEQ.f CHEMEV.f CHMIMP.f CHPACK.f \
	DIFFx.f DIFFY.f DIFFZ.f DUMP.f EXITA.f FRAN.f HOT.f INJECT.f \
	LAVA.f LNPKD.f MOMFLX.f MOMFLY.f MOMFLZ.f MOMSRC.f MOMX.f MOMY.f \
	MOMZ.f OLDTIM.f PDVETC.f PFIND.f PINT.f PLOT.f TECPLOT.f HDFOUT.f \
        PMOVE.f PARTIHEAT.f\
	PRINT.f PSTATE.f RDIVID.f REPACK.f RINPUT.f RMULT.f RSTART.f \
	SETQ.f SETR.f SETUP.f START.f STATE.f TIMSTP.f USER.f \
	USERPP.f WALLF.f VISCR.f VISCX.f VISCY.f VISCZ.f XPORT.f

HEADERS = COML.h COMMPTC.h

#LIBS = -lvec
LIBS =

LFLAGS =

INCS =

FC = f77

FFLAGS = -O -u -w 
#FFLAGS = -O -u -w -pg
#FFLAGS = -O +OPP -u -w
#FFLAGS = -O +OP -Wl,-a,archive -u -w

#PROG = ../data/LAVA_Metco_9MB
PROG = lava.x

.SUFFIXES:

.SUFFIXES: .f .o

.f.o:
#	${FC} ${FFLAGS} ${INCS} -WP,-l=`basename $*.f .f`.l,-listoptions=o -c $*.f
	${FC} ${FFLAGS} ${INCS} -c $*.f

${PROG}: ${SRCS:.f=.o}
	${FC} ${FFLAGS} ${INCS} -o $@ ${SRCS:.f=.o} ${LFLAGS} ${LIBS}

clean:
	rm -f ${SRCS:.f=.o} ${PROG} core

depend:
	makedepend -- ${HEADERS} -- ${SRCS}

# DO NOT DELETE THIS LINE -- make depend depends on it.

${SRCS:.f=.o} : COML.h COMMPTC.h
