#pragma once
#include "afxpropertygridctrl.h"
class CMyValueProperty :
	public CMFCPropertyGridProperty
{
	friend class CMFCPropertyGridCtrl;

	DECLARE_DYNAMIC(CMyValueProperty)

	 //Construction
public:
	CMyValueProperty(const CString& strName, UINT varValue, UINT maxValue, UINT minValue, LPCTSTR lpszDescr = NULL) :
		CMFCPropertyGridProperty(strName, (SHORT)varValue, lpszDescr)
	{
		m_max_value = maxValue;
		m_min_value = minValue;
	}
	CMyValueProperty(const CString& strName, CString varValue, LPCTSTR lpszDescr = NULL) :
		CMFCPropertyGridProperty(strName, varValue, lpszDescr)
	{
		
	}

	virtual ~CMyValueProperty()
	{}

	virtual BOOL OnEndEdit();
	// Attributes
public:

	// Attributes
protected:
	UINT m_max_value;
	UINT m_min_value;
};

