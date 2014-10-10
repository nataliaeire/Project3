TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    includingClasses.cpp \
    initialTwoBodySystem.cpp \
    main.cpp \
    rk4.cpp \
    system.cpp \
    celestialbody.cpp \
    vec3.cpp
LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS += -Wall

HEADERS += \
    rk4.h \
    system.h \
    celestialbody.h \
    vec3.h
