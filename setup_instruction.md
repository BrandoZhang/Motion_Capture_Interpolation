Code setup instructions
====================================================================================
How to build fltk in Windows, Linux, MAC OS X
====================================================================================
Please follow the corresponding readme files in the folder: fltk-1.3.8.
Note: On Linux, you may need to install the "autoconf" tool.

Note: There are many projects in FLTK library. This homework
is only dependent on fltk and fltkgl projects. Students can choose
to only compile those two projects in the Visual Studio project
of the FLTK library.

Note: If you want to compile the homework under the release
configuration, you should first compile the FLTK projects in release.
Same goes for the debug configuration.

Note: Both the main homework project and FLTK projects work only
on the Win32 platform. Please do not use other platforms (e.g. X64).

====================================================================================
How to build starter code in Linux, MAC OS X (Assuming fltk has been compiled)
====================================================================================
1) Enter the %homeworkFolder%/mocapPlayer-starter
2) make

====================================================================================
How to build starter code using Visual Studio 2017 (Assuming fltk has been compiled)
====================================================================================
1) Open the project file in homework folder: IDE-starter/VS2017/mocapPlayer.sln
2) Choose Debug/Release mode
3) Compile project: mocapPlayer
4) Compile project: interpolate

Because there are many versions of Windows 10 (different build versions),
you may get an error when compiling: ... selecting "Retarget solution".
This can be solved, as the message says, by right-clicking the solution
and selecting "Retarget solution".