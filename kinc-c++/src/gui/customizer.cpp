#include "customizer.h"
#include <QObject>
#include <QAction>
#include <QMenu>
#include "../core/datafactory.h"
#include "../core/analyticfactory.h"



void Customizer::setupDataMenu(QMenu &menu, QObject *parent, const std::function<void(int)>& callback)
{
    EAbstractDataFactory& factory {DataFactory::instance()};
    for (quint16 i = 0; i < factory.size() ;++i)
    {
        QAction* action {new QAction(factory.name(i),parent)};
        QObject::connect(action,&QAction::triggered,[callback,i]{ callback(i); });
        menu.addAction(action);
    }
}



void Customizer::setupAnalyticMenu(QMenu &menu, QObject *parent, const std::function<void(int)>& callback)
{
    EAbstractAnalyticFactory& factory {AnalyticFactory::instance()};
    for (quint16 i = 0; i < factory.size() ;++i)
    {
        QAction* action {new QAction(factory.name(i),parent)};
        QObject::connect(action,&QAction::triggered,[callback,i]{ callback(i); });
        menu.addAction(action);
    }
}



QString Customizer::aboutTitle()
{
    return QStringLiteral("About KINC");
}



QString Customizer::aboutRichText()
{
    return QStringLiteral("<p>Knowledge Independent Network Construction (KINC) is a high-performance application for constructing Gene Co-expression Networks (GCNs) from gene expression data.</p>");
}



QPixmap Customizer::splashImage()
{
    return QPixmap();
}



QIcon Customizer::icon()
{
    return QIcon();
}
