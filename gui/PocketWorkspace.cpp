#include "PocketWorkspace.h"
#include "ligint/pocket_interfaces.hpp"
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QListWidget>
#include <QLabel>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QGroupBox>
#include <QPainter>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QtConcurrent/QtConcurrent>
#include <algorithm>
#include <cmath>

namespace gui {
class PocketCanvas : public QWidget {
public:
    explicit PocketCanvas(PocketWorkspace* owner):QWidget(owner),owner_(owner){
        setMinimumSize(400,300);setToolTip("Drag to rotate; wheel to zoom");
    }
protected:
    void mousePressEvent(QMouseEvent* e) override {last_=e->pos();}
    void mouseMoveEvent(QMouseEvent* e) override {
        if(e->buttons()&Qt::LeftButton){yaw_+=(e->pos().x()-last_.x())*0.008;
            pitch_+=(e->pos().y()-last_.y())*0.008;update();}
        last_=e->pos();
    }
    void wheelEvent(QWheelEvent* e) override {
        zoom_=std::clamp(zoom_*std::pow(1.001,e->angleDelta().y()),0.2,10.0);update();
    }
    void paintEvent(QPaintEvent*) override {
        QPainter p(this);p.fillRect(rect(),QColor(16,23,34));p.setRenderHint(QPainter::Antialiasing);
        const auto members=owner_->visibleMembers();
        if(members.empty()){p.setPen(Qt::white);p.drawText(rect(),Qt::AlignCenter,"Run Pocket Finder and select a candidate");return;}
        std::vector<ligint::PocketPoint> points;
        ligint::V3 center{0,0,0};
        for(auto id:members)if(id<owner_->result_.pockets.size()){
            const auto& pocket=owner_->result_.pockets[id];
            points.insert(points.end(),pocket.points.begin(),pocket.points.end());
        }
        if(points.empty())return;
        for(const auto& a:points)for(int k=0;k<3;k++)center[k]+=a.pos[k]/points.size();
        double extent=1;
        for(const auto& a:points)for(int k=0;k<3;k++)extent=std::max(extent,std::fabs(a.pos[k]-center[k]));
        const double scale=0.43*std::min(width(),height())*zoom_/extent;
        struct Dot {double depth;QPointF screen;QColor color;};
        std::vector<Dot> dots;dots.reserve(points.size());
        static const QColor palette[]={QColor(40,185,230),QColor(255,160,60),QColor(120,220,130),QColor(214,100,236)};
        const double cy=std::cos(yaw_),sy=std::sin(yaw_),cp=std::cos(pitch_),sp=std::sin(pitch_);
        for(std::size_t m=0;m<members.size();m++){
            const auto id=members[m];if(id>=owner_->result_.pockets.size())continue;
            for(const auto& voxel:owner_->result_.pockets[id].points){
                const double x=voxel.pos[0]-center[0],y=voxel.pos[1]-center[1],z=voxel.pos[2]-center[2];
                const double rx=cy*x+sy*z,rz=-sy*x+cy*z;
                const double ry=cp*y-sp*rz,depth=sp*y+cp*rz;
                QColor color=palette[m%4];
                color.setAlphaF(std::clamp(0.18+0.65*voxel.buriedness,0.0,1.0));
                dots.push_back({depth,QPointF(width()/2.0+rx*scale,height()/2.0-ry*scale),color});
            }
        }
        std::sort(dots.begin(),dots.end(),[](const Dot&a,const Dot&b){return a.depth<b.depth;});
        p.setPen(Qt::NoPen);
        for(const auto& d:dots){p.setBrush(d.color);p.drawEllipse(d.screen,2.5,2.5);}
        p.setPen(QColor(230,238,248));
        p.drawText(QRect(10,8,width()-20,24),Qt::AlignLeft,
            QString("Pocket voxels: %1  |  drag: rotate  wheel: zoom").arg(points.size()));
    }
private:
    PocketWorkspace* owner_;QPoint last_;double yaw_=0.5,pitch_=0.4,zoom_=1;
};

PocketWorkspace::PocketWorkspace(QWidget* parent):QWidget(parent){
    auto* outer=new QHBoxLayout(this);
    auto* controls=new QWidget;
    controls->setMaximumWidth(290);
    auto* column=new QVBoxLayout(controls);
    auto* form=new QFormLayout;
    spacing_=new QDoubleSpinBox;spacing_->setRange(0.5,3);spacing_->setValue(1);spacing_->setSingleStep(.25);
    buried_=new QDoubleSpinBox;buried_->setRange(.1,1);buried_->setSingleStep(.05);buried_->setValue(.5);
    ray_=new QDoubleSpinBox;ray_->setRange(2,30);ray_->setValue(10);
    surface_=new QDoubleSpinBox;surface_->setRange(1,15);surface_->setValue(5);
    probe_=new QDoubleSpinBox;probe_->setRange(.5,3);probe_->setValue(1.4);
    radius_=new QSpinBox;radius_->setRange(0,2);radius_->setValue(1);
    minPts_=new QSpinBox;minPts_->setRange(1,500);minPts_->setValue(30);
    morph_=new QComboBox;morph_->addItems({"none","opening","closing","open-close"});morph_->setCurrentText("opening");
    form->addRow("Grid spacing (Å)",spacing_);form->addRow("Min buriedness",buried_);
    form->addRow("Ray length (Å)",ray_);form->addRow("Surface distance (Å)",surface_);
    form->addRow("Solvent probe (Å)",probe_);form->addRow("Morphology",morph_);
    form->addRow("Morph radius",radius_);form->addRow("Min cluster voxels",minPts_);
    column->addLayout(form);
    run_=new QPushButton("Find pockets");column->addWidget(run_);
    status_=new QLabel("Open a structure to start");status_->setWordWrap(true);column->addWidget(status_);
    display_=new QComboBox;display_->addItems({"Subpockets","Candidate groups"});
    column->addWidget(display_);
    entries_=new QListWidget;column->addWidget(entries_,1);
    auto* canvas=new PocketCanvas(this);canvas_=canvas;
    outer->addWidget(controls);outer->addWidget(canvas,1);
    connect(display_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this]{
        entries_->clear();
        if(display_->currentIndex()==0){
            for(std::size_t i=0;i<result_.pockets.size();++i)
                entries_->addItem(QString("#%1 | %2 Å³ | buried %3").arg(i+1)
                     .arg(result_.pockets[i].volume,0,'f',0)
                     .arg(result_.pockets[i].meanBuriedness,0,'f',2));
        }else{
            for(const auto& g:groups_){
                QStringList ids;for(auto id:g.members)ids<<QString::number(id+1);
                entries_->addItem(QString("[%1] | %2 Å³ | score %3")
                    .arg(ids.join("+")).arg(g.volume,0,'f',0).arg(g.score,0,'f',3));
            }
        }
        if(entries_->count())entries_->setCurrentRow(0);
        canvas_->update();
    });
    connect(entries_,&QListWidget::currentRowChanged,canvas_,QOverload<>::of(&QWidget::update));
    connect(run_,&QPushButton::clicked,this,[this]{
        if(path_.isEmpty())return;
        ligint::PocketOptions o;
        o.spacing=spacing_->value();o.minBuriedness=buried_->value();o.rayLength=ray_->value();
        o.maxSurfaceDistance=surface_->value();o.solventProbeRadius=probe_->value();
        o.morphology=morph_->currentText().toStdString();o.morphologyRadius=radius_->value();
        o.minClusterPoints=static_cast<std::size_t>(minPts_->value());
        const QString path=path_;
        run_->setEnabled(false);status_->setText("Computing pockets...");
        watcher_.setFuture(QtConcurrent::run([path,o]{
            PocketGuiOutcome out;
            try{
                auto atoms=ligint::readPocketAtoms(path.toStdString());
                out.result=ligint::findPockets(atoms,o);
                std::stable_sort(out.result.pockets.begin(),out.result.pockets.end(),
                    [](const auto&a,const auto&b){return a.heuristicScore>b.heuristicScore;});
                const auto edges=ligint::analyzePocketInterfaces(out.result.pockets,atoms,o.spacing);
                out.groups=ligint::rankPocketGroups(out.result.pockets,edges);
            }catch(const std::exception& e){out.error=QString::fromUtf8(e.what());}
            return out;
        }));
    });
    connect(&watcher_,&QFutureWatcher<PocketGuiOutcome>::finished,this,[this]{
        auto out=watcher_.result();run_->setEnabled(true);
        if(!out.error.isEmpty()){status_->setText("Error: "+out.error);return;}
        result_=std::move(out.result);groups_=std::move(out.groups);
        status_->setText(QString("%1 pockets, %2 groups; %3 atoms")
            .arg(result_.pockets.size()).arg(groups_.size()).arg(result_.selectedAtoms));
        display_->setCurrentIndex(1);display_->setCurrentIndex(0);
        canvas_->update();
    });
}
PocketWorkspace::~PocketWorkspace(){if(watcher_.isRunning())watcher_.waitForFinished();}
void PocketWorkspace::setStructurePath(const QString& path){
    path_=path;run_->setEnabled(!path.isEmpty());
    status_->setText(path.isEmpty()?"Open a structure to start":"Ready: "+path);
    result_={};groups_.clear();entries_->clear();canvas_->update();
}
std::vector<std::size_t> PocketWorkspace::visibleMembers() const{
    const int row=entries_->currentRow();
    if(row<0)return {};
    if(display_->currentIndex()==0){
        return static_cast<std::size_t>(row)<result_.pockets.size()
            ?std::vector<std::size_t>{static_cast<std::size_t>(row)}:std::vector<std::size_t>{};
    }
    return static_cast<std::size_t>(row)<groups_.size()?groups_[row].members:std::vector<std::size_t>{};
}
} // namespace gui
