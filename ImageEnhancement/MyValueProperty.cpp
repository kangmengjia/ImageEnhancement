#include "stdafx.h"
#include "MyValueProperty.h"
#include "MainFrm.h"
#include "resource.h"
#include "ImageEnhancementView.h"

IMPLEMENT_DYNAMIC(CMyValueProperty, CMFCPropertyGridProperty)

BOOL CMyValueProperty::OnEndEdit(){
	if (m_varValue.vt == VT_I2){
		if (m_varValue.iVal < m_min_value){
			m_varValue.iVal = (SHORT)m_min_value;
		}
		else if (m_varValue.iVal > m_max_value){
			m_varValue.iVal = (SHORT)m_max_value;
		}
	}
	((CImageEnhancementView*)(((CMainFrame*)AfxGetMainWnd())->GetActiveView()))->PostMessageW(ID_RELOAD_MSG, MAKEWPARAM(ID_RELOAD_MSG, BN_CLICKED), (LPARAM)this->GetParent());
	return CMFCPropertyGridProperty::OnEndEdit();
}