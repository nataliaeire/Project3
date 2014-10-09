TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    includingClasses.cpp \
    initialTwoBodySystem.cpp \
    KeplerClasses.cpp \
    main.cpp
LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS += -Wall

HEADERS += \
    KeplerClasses.h
