function caseSettings = getSystemSettings(sysName,caseName)
caseSettings = eval([sysName '_settings(caseName)']);
end