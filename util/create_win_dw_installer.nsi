; Requires https://nsis.sourceforge.io/EnVar_plug-in
; just extract to C:\Program (x86)\NSIS

Name "deconwolf"

;Icon "todo.ico"

OutFile "deconwolf_0_4_0_Setup.exe"

; The default installation directory
InstallDir $PROGRAMFILES\deconwolf

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\deconwolf" "Install_Dir"

RequestExecutionLevel admin

!include LogicLib.nsh

;--------------------------------

; Pages

Page components
Page directory
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

;--------------------------------

; The stuff to install
Section "deconwolf"
 
  SectionIn RO
 
   ; Add the istall dir to the path
  EnVar::SetHKCU
EnVar::Check "Path" "$InstDir"
Pop $0
${If} $0 = 0
  DetailPrint "Already there"
${Else}
  EnVar::AddValue "Path" "$InstDir"
  Pop $0 ; 0 on success
${EndIf}
 
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  
  ; Put file there (you can add more File lines too)
  File "dw.exe"
  File "dw_bw.exe"
  File "fftw3f.dll"
  File "getopt.dll"
  File "gsl.dll"
  File "gslcblas.dll"
  File "jpeg62.dll"
  File "liblzma.dll"
  File "libomp.dll"
  File "tiff.dll"
  File "zlib1.dll"

  ; NOTE: Same list for uninstaller
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\deconwolf "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\deconwolf" "DisplayName" "deconwolf"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\deconwolf" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\deconwolf" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\deconwolf" "NoRepair" 1
  WriteUninstaller "$INSTDIR\uninstall.exe"
  
SectionEnd

; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts (required)"
  SectionIn RO

  CreateDirectory "$SMPROGRAMS\deconwolf"
  CreateShortcut "$SMPROGRAMS\deconwolf\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0  
  
SectionEnd

;--------------------------------

; Uninstaller

Section "Uninstall"
  
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\deconwolf"
  DeleteRegKey HKLM SOFTWARE\deconwolf
  
  Delete "$INSTDIR\uninstall.exe"
  Delete "$INSTDIR\dw.exe"
  Delete "$INSTDIR\dw_bw.exe"
  Delete "$INSTDIR\fftw3f.dll"
  Delete "$INSTDIR\getopt.dll"
  Delete "$INSTDIR\gsl.dll"
  Delete "$INSTDIR\gslcblas.dll"
  Delete "$INSTDIR\jpeg62.dll"
  Delete "$INSTDIR\liblzma.dll"
  Delete "$INSTDIR\libomp.dll"
  Delete "$INSTDIR\tiff.dll"
  Delete "$INSTDIR\zlib1.dll"

  RMDir "$INSTDIR"

  ; Remove from path
  EnVar::DeleteValue "Path" "$InstDir"
Pop $0

SectionEnd
