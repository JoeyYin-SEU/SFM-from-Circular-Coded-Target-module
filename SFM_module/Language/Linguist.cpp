
#include "Linguist.h"
#include <qcoreapplication.h>
CLinguist* CLinguist::linguist = nullptr;
CLinguist::CLinguist()
{
    m_trans = new QTranslator;
}

CLinguist* CLinguist::GetLinguistInstance()
{
    if (nullptr == linguist)
    {
        linguist = new CLinguist;
    }

    return linguist;
}

void CLinguist::ChangeLanguage(Language lan)
{
    if (lan == this->m_CurrentLanguage) return;
    bool ret;
    switch (lan)
    {
    case zh_cn:
        if (nullptr != m_trans)
        {
            qApp->removeTranslator(m_trans);
        }
        ret = m_trans->load("pet_zh.qm");
        if (ret)
        {
            qApp->installTranslator(m_trans);
        }
        break;
    case en_us:
        if (nullptr != m_trans)
        {
            qApp->removeTranslator(m_trans);
        }
        ret = m_trans->load("pet_en.qm");
        if (ret)
        {
            qApp->installTranslator(m_trans);
        }
        break;
    case zhc_cn:
        if (nullptr != m_trans)
        {
            qApp->removeTranslator(m_trans);
        }
        ret = m_trans->load("pet_zhc.qm");
        if (ret)
        {
            qApp->installTranslator(m_trans);
        }
        break;
    default:
        break;
    }
    //if (ret)
    {
        this->m_CurrentLanguage = lan;
        // 发出语言被切换的信号
        emit LanguageChaned();
    }
}