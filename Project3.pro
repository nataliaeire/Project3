TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    includingClasses.cpp \
    initialTwoBodySystem.cpp \
    main.cpp \
    system.cpp \
    celestialbody.cpp \
    vec3.cpp \
    printing.cpp \
    integrator.cpp \
    gaussiandeviate.cpp \
    cpelapsedtimer.cpp \
    oldmain.cpp
#LIBS += -larmadillo -llapack -lblas

#QMAKE_CXXFLAGS += -Wall

HEADERS += \
    system.h \
    celestialbody.h \
    vec3.h \
    integrator.h \
    printing.h \
    gaussiandeviate.h \
    cpelapsedtimer.h

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS   += -fopenmp
