// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef GLOBALATTRIBUTES_H
#define GLOBALATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: GlobalAttributes
//
// Purpose:
//    This class contains attributes associated with the main window.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//
// ****************************************************************************

class STATE_API GlobalAttributes : public AttributeSubject
{
public:
    enum PrecisionType
    {
        Float,
        Native,
        Double
    };
    enum BackendType
    {
        VTK,
        VTKM
    };

    // These constructors are for objects of this class
    GlobalAttributes();
    GlobalAttributes(const GlobalAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    GlobalAttributes(private_tmfs_t tmfs);
    GlobalAttributes(const GlobalAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~GlobalAttributes();

    virtual GlobalAttributes& operator = (const GlobalAttributes &obj);
    virtual bool operator == (const GlobalAttributes &obj) const;
    virtual bool operator != (const GlobalAttributes &obj) const;
private:
    void Init();
    void Copy(const GlobalAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSources();
    void SelectWindows();

    // Property setting methods
    void SetSources(const stringVector &sources_);
    void SetWindows(const intVector &windows_);
    void SetActiveWindow(int activeWindow_);
    void SetIconifiedFlag(bool iconifiedFlag_);
    void SetAutoUpdateFlag(bool autoUpdateFlag_);
    void SetReplacePlots(bool replacePlots_);
    void SetApplyOperator(bool applyOperator_);
    void SetApplySelection(bool applySelection_);
    void SetApplyWindow(bool applyWindow_);
    void SetExecuting(bool executing_);
    void SetWindowLayout(int windowLayout_);
    void SetMakeDefaultConfirm(bool makeDefaultConfirm_);
    void SetCloneWindowOnFirstRef(bool cloneWindowOnFirstRef_);
    void SetAutomaticallyAddOperator(bool automaticallyAddOperator_);
    void SetTryHarderCyclesTimes(bool tryHarderCyclesTimes_);
    void SetTreatAllDBsAsTimeVarying(bool treatAllDBsAsTimeVarying_);
    void SetCreateMeshQualityExpressions(bool createMeshQualityExpressions_);
    void SetCreateTimeDerivativeExpressions(bool createTimeDerivativeExpressions_);
    void SetCreateVectorMagnitudeExpressions(bool createVectorMagnitudeExpressions_);
    void SetNewPlotsInheritSILRestriction(bool newPlotsInheritSILRestriction_);
    void SetUserDirForSessionFiles(bool userDirForSessionFiles_);
    void SetSaveCrashRecoveryFile(bool saveCrashRecoveryFile_);
    void SetIgnoreExtentsFromDbs(bool ignoreExtentsFromDbs_);
    void SetExpandNewPlots(bool expandNewPlots_);
    void SetUserRestoreSessionFile(bool userRestoreSessionFile_);
    void SetPrecisionType(PrecisionType precisionType_);
    void SetBackendType(BackendType backendType_);
    void SetRemoveDuplicateNodes(bool removeDuplicateNodes_);

    // Property getting methods
    const stringVector &GetSources() const;
          stringVector &GetSources();
    const intVector    &GetWindows() const;
          intVector    &GetWindows();
    int                GetActiveWindow() const;
    bool               GetIconifiedFlag() const;
    bool               GetAutoUpdateFlag() const;
    bool               GetReplacePlots() const;
    bool               GetApplyOperator() const;
    bool               GetApplySelection() const;
    bool               GetApplyWindow() const;
    bool               GetExecuting() const;
    int                GetWindowLayout() const;
    bool               GetMakeDefaultConfirm() const;
    bool               GetCloneWindowOnFirstRef() const;
    bool               GetAutomaticallyAddOperator() const;
    bool               GetTryHarderCyclesTimes() const;
    bool               GetTreatAllDBsAsTimeVarying() const;
    bool               GetCreateMeshQualityExpressions() const;
    bool               GetCreateTimeDerivativeExpressions() const;
    bool               GetCreateVectorMagnitudeExpressions() const;
    bool               GetNewPlotsInheritSILRestriction() const;
    bool               GetUserDirForSessionFiles() const;
    bool               GetSaveCrashRecoveryFile() const;
    bool               GetIgnoreExtentsFromDbs() const;
    bool               GetExpandNewPlots() const;
    bool               GetUserRestoreSessionFile() const;
    PrecisionType      GetPrecisionType() const;
    BackendType        GetBackendType() const;
    bool               GetRemoveDuplicateNodes() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string PrecisionType_ToString(PrecisionType);
    static bool PrecisionType_FromString(const std::string &, PrecisionType &);
protected:
    static std::string PrecisionType_ToString(int);
public:
    static std::string BackendType_ToString(BackendType);
    static bool BackendType_FromString(const std::string &, BackendType &);
protected:
    static std::string BackendType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_sources = 0,
        ID_windows,
        ID_activeWindow,
        ID_iconifiedFlag,
        ID_autoUpdateFlag,
        ID_replacePlots,
        ID_applyOperator,
        ID_applySelection,
        ID_applyWindow,
        ID_executing,
        ID_windowLayout,
        ID_makeDefaultConfirm,
        ID_cloneWindowOnFirstRef,
        ID_automaticallyAddOperator,
        ID_tryHarderCyclesTimes,
        ID_treatAllDBsAsTimeVarying,
        ID_createMeshQualityExpressions,
        ID_createTimeDerivativeExpressions,
        ID_createVectorMagnitudeExpressions,
        ID_newPlotsInheritSILRestriction,
        ID_userDirForSessionFiles,
        ID_saveCrashRecoveryFile,
        ID_ignoreExtentsFromDbs,
        ID_expandNewPlots,
        ID_userRestoreSessionFile,
        ID_precisionType,
        ID_backendType,
        ID_removeDuplicateNodes,
        ID__LAST
    };

private:
    stringVector sources;
    intVector    windows;
    int          activeWindow;
    bool         iconifiedFlag;
    bool         autoUpdateFlag;
    bool         replacePlots;
    bool         applyOperator;
    bool         applySelection;
    bool         applyWindow;
    bool         executing;
    int          windowLayout;
    bool         makeDefaultConfirm;
    bool         cloneWindowOnFirstRef;
    bool         automaticallyAddOperator;
    bool         tryHarderCyclesTimes;
    bool         treatAllDBsAsTimeVarying;
    bool         createMeshQualityExpressions;
    bool         createTimeDerivativeExpressions;
    bool         createVectorMagnitudeExpressions;
    bool         newPlotsInheritSILRestriction;
    bool         userDirForSessionFiles;
    bool         saveCrashRecoveryFile;
    bool         ignoreExtentsFromDbs;
    bool         expandNewPlots;
    bool         userRestoreSessionFile;
    int          precisionType;
    int          backendType;
    bool         removeDuplicateNodes;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define GLOBALATTRIBUTES_TMFS "s*i*ibbbbbbbibbbbbbbbbbbbbbiib"

#endif