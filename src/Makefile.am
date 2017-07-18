#
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS = oemtp
oemtp_SOURCES = oemtp.cpp

lib_LTLIBRARIES = liboemtp-@OEMTP_API_VERSION@.la
liboemtp_@OEMTP_API_VERSION@_la_SOURCES = analysis.cpp\
components.cpp\
cable.cpp\
interface_orelay.cpp\
ode2.cpp\
block.cpp\
transformers.cpp\
support/support.cpp\
support/fdump.f\
support/d1mach.f\
support/xerbla.f\
support/zacon.f\
support/zdiv.f\
support/zlacgv.f\
support/zlaqr3.f\
support/zlaset.f\
support/zseri.f\
support/zunhj.f\
support/zuoik.f\
support/dgamln.f\
support/i1mach.f\
support/xercnt.f\
support/zairy.f\
support/zexp.f\
support/zlacpy.f\
support/zlaqr4.f\
support/zlassq.f\
support/zshch.f\
support/zuni1.f\
support/zwrsk.f\
support/disnan.f\
support/ieeeck.f\
support/xerhlt.f\
support/zasyi.f\
support/zgebak.f\
support/zladiv.f\
support/zlaqr5.f\
support/zlatrs.f\
support/zsqrt.f\
support/zuni2.f\
support/dlabad.f\
support/ilaenv.f\
support/xermsg.f\
support/zbesi.f\
support/zgebal.f\
support/zlahqr.f\
support/zlarfb.f\
support/zlog.f\
support/ztrevc.f\
support/zunik.f\
support/dladiv.f\
support/ilazlc.f\
support/xerprn.f\
support/zbesk.f\
support/zgeev.f\
support/zlahr2.f\
support/zlarf.f\
support/zmlri.f\
support/ztrexc.f\
support/zunk1.f\
support/dlaisnan.f\
support/ilazlr.f\
support/xersve.f\
support/zbinu.f\
support/zgehd2.f\
support/zlange.f\
support/zlarfg.f\
support/zmlt.f\
support/zuchk.f\
support/zunk2.f\
support/dlamch.f\
support/iparmq.f\
support/xgetua.f\
support/zbknu.f\
support/zgehrd.f\
support/zlaqr0.f\
support/zlarft.f\
support/zrati.f\
support/zung2r.f\
support/zunm2r.f\
support/dlapy2.f\
support/j4save.f\
support/zabs.f\
support/zbuni.f\
support/zhseqr.f\
support/zlaqr1.f\
support/zlartg.f\
support/zrot.f \
support/zunghr.f\
support/zunmhr.f\
support/dlapy3.f\
support/lsame.f\
support/zacai.f\
support/zbunk.f\
support/zkscl.f\
support/zlaqr2.f\
support/zlascl.f\
support/zs1s2.f\
support/zungqr.f\
support/zunmqr.f
#components.cpp \
#interface_orelay.cpp


liboemtp_@OEMTP_API_VERSION@_la_LDFLAGS = -version-info $(OEMTP_SO_VERSION)
oemtp_includedir = $(includedir)/oemtp-$(OEMTP_API_VERSION)
nobase_oemtp_include_HEADERS = analysis.h\
components.h\
transformers.h\
cable.h\
constants.h\
interface_orelay.h\
ode2.h\
block.h\
support/support.h
#components.h \
#interface_orelay.h

oemtp_libincludedir = $(libdir)/oemtp-$(OEMTP_API_VERSION)/include

oemtp_LDADD = ./liboemtp-${OEMTP_API_VERSION}.la










