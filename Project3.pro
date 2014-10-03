TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    includingClasses.cpp
LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS += -Wall
