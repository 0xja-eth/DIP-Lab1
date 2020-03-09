#include "stdafx.h"
#include "CEditEx.h"

// CEditEx

IMPLEMENT_DYNAMIC(CEditEx, CEdit)

CEditEx::CEditEx() {}

CEditEx::~CEditEx() {}

BEGIN_MESSAGE_MAP(CEditEx, CEdit)
	ON_WM_CHAR()
	ON_WM_KEYUP()
END_MESSAGE_MAP()

// CEditEx 消息处理程序
void CEditEx::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags) {
	// TODO:  在此添加消息处理程序代码和/或调用默认值

	// 处理小数点
	if (nChar == '.') {
		CString str;
		GetWindowText(str);

		// 限制第一位为小数
		if (str.GetLength() == 0) {
			// 第一位输入小数点
			MessageBox(_T("第一位不可以是小数点"));
			return;
		}
		// 限制只允许有一个小数点
		if (str.Find('.') == -1) {
			CEdit::OnChar(nChar, nRepCnt, nFlags);
		} else {
			CString str;
			GetWindowText(str);
			if (str[0] == '.') {
				SetWindowText(str.Mid(1, str.GetLength()));
				MessageBox(_T("第一位不可以是小数点"));
			}
			// 小数点出现第二次
			MessageBox(_T("小数点只能输入一次"));
		}
	}
	// 数理数字和退格键
	else if ((nChar >= '0' && nChar <= '9') || nChar == 0x08) {
		CEdit::OnChar(nChar, nRepCnt, nFlags);
	} else {
		// 出现非数字键，退格键
		MessageBox(_T("只能输入数字，退格键"));
	}
}

// 修复先输入数字之后，可以在第一位输入小数点
void CEditEx::OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags) {
	// TODO: 在此添加消息处理程序代码和/或调用默认值
	if (nChar == VK_DECIMAL || nChar == VK_OEM_PERIOD) {
		CString str;
		GetWindowText(str);
		if (str[0] == '.') {
			SetWindowText(str.Mid(1, str.GetLength()));
			MessageBox(_T("第一位不可以是小数点"));
		}
	}
	CEdit::OnKeyUp(nChar, nRepCnt, nFlags);

}
