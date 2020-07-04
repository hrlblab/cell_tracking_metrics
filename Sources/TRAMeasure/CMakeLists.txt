########################
#  Application settings
########################

cmake_minimum_required(VERSION 3.3)

########################
#  Compilation flags
########################

IF (WIN32)

	ADD_DEFINITIONS (-D_CRT_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_DEPRECATE -DWXUSINGDLL)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /nologo /wd4251 /wd4275")

ENDIF (WIN32)

##########################################################################
# find required libraries
##########################################################################

# wxWidgets
########################

IF (WIN32)

	FIND_LIBRARY (LIB_WXBASE wxbase31u DOC "wxbase31u.lib file")		
	SET(WX_LIBS ${LIB_WXBASE})

	# set path to setup.h - the configuration file of wxWidgets libraries
	GET_FILENAME_COMPONENT (LIB_WX_DIR ${LIB_WXBASE} PATH)	
	INCLUDE_DIRECTORIES (${LIB_WX_DIR}/mswu)

	# set path to headers of wxWidgets
	FIND_PATH (WX_HEADERS wx/wx.h ${INC_DIRS} DOC "Path to wxWidgets header files")
	INCLUDE_DIRECTORIES (${WX_HEADERS})

ENDIF (WIN32)


# i3dlibs
########################

FIND_LIBRARY (LIB_I3DCORE i3dcore DOC "i3dcore.lib file")
FIND_LIBRARY (LIB_I3DALGO i3dalgo DOC "i3dalgo.lib file")

SET(I3D_LIBS ${LIB_I3DCORE} ${LIB_I3DALGO})

# set path to headers of i3dlibs
FIND_PATH (I3D_HEADERS i3d/image3d.h ${INC_DIRS} DOC "Path to i3dlibs header files")
INCLUDE_DIRECTORIES (${I3D_HEADERS})

# build-up the application
ADD_EXECUTABLE(TRAMeasure TRAMeasure.cpp)
TARGET_LINK_LIBRARIES (TRAMeasure ${I3D_LIBS} ${WX_LIBS})
