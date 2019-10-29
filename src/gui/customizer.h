#ifndef CUSTOMIZER_H
#define CUSTOMIZER_H
#include <ace/gui/eabstractcustomizer.h>



class Customizer : public EAbstractCustomizer
{
    virtual void setupDataMenu(QMenu& menu, QObject* parent, const std::function<void(int)>& callback) override final;
    virtual void setupAnalyticMenu(QMenu& menu, QObject* parent, const std::function<void(int)>& callback) override final;
    virtual QString aboutTitle() override final;
    virtual QString aboutRichText() override final;
    virtual QPixmap splashImage() override final;
    virtual QIcon icon() override final;
};



#endif
