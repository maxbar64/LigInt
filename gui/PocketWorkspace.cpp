#include "PocketWorkspace.h"
#include "ligint/pocket_interfaces.hpp"
#include "PocketSmoothMesh.h"
#include "CommandLine.h"
#include <QComboBox>
#include <QCheckBox>
#include <QApplication>
#include <QClipboard>
#include <QStringList>
#include <QPen>
#include <QAbstractItemView>
#include <QSignalBlocker>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include <QPolygonF>
#include <QFrame>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QListWidget>
#include <QLabel>
#include <QPushButton>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QGroupBox>
#include <QScrollArea>
#include <QStackedWidget>
#include <cstdio>
#include <iomanip>
#include <QToolButton>
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
        if(owner_->algorithm_->currentIndex()==1 && owner_->alphaRender_->currentIndex()==1) {
            if(members.empty()){p.setPen(Qt::white);p.drawText(rect(),Qt::AlignCenter,"Run Alpha-sphere search and select a candidate");return;}
            struct AlphaDraw {const ligint::AlphaSphere* sphere;std::size_t pocket;};
            std::vector<AlphaDraw> spheres;
            for(const auto id:members)if(id<owner_->alphaPockets_.size())
                for(const auto& sphere:owner_->alphaPockets_[id].spheres)spheres.push_back({&sphere,id});
            if(spheres.empty())return;
            ligint::V3 center{0,0,0};
            for(const auto& item:spheres)for(int k=0;k<3;++k)center[k]+=item.sphere->center[k]/spheres.size();
            double extent=1;
            for(const auto& item:spheres)for(int k=0;k<3;++k)
                extent=std::max(extent,std::fabs(item.sphere->center[k]-center[k])+item.sphere->radius);
            const double scale=0.43*std::min(width(),height())*zoom_/extent;
            const double cy=std::cos(yaw_),sy=std::sin(yaw_),cp=std::cos(pitch_),sp=std::sin(pitch_);
            auto project=[&](const ligint::V3& v){
                const double x=v[0]-center[0],y=v[1]-center[1],z=v[2]-center[2];
                const double rx=cy*x+sy*z,rz=-sy*x+cy*z;
                return std::pair<QPointF,double>{QPointF(width()/2.0+rx*scale,height()/2.0-(cp*y-sp*rz)*scale),sp*y+cp*rz};
            };
            for(std::size_t i=0;i<owner_->atoms_.size();i+=2){
                const auto& atom=owner_->atoms_[i];double d2=0;
                for(int k=0;k<3;++k)d2+=std::pow(atom.pos[k]-center[k],2);
                if(d2>400)continue;
                const auto screen=project(atom.pos);
                p.setPen(Qt::NoPen);p.setBrush(QColor(148,159,171,75));p.drawEllipse(screen.first,1.3,1.3);
            }
            static const QColor colors[]={QColor(40,185,230),QColor(255,160,60),QColor(120,220,130),QColor(214,100,236)};
            std::sort(spheres.begin(),spheres.end(),[&](const auto& a,const auto& b){
                return project(a.sphere->center).second<project(b.sphere->center).second;
            });
            for(const auto& item:spheres){
                const auto screen=project(item.sphere->center);
                QColor fill=colors[item.pocket%4];fill.setAlpha(std::clamp(owner_->alphaOpacity_->value()/3,20,90));
                QColor edge=colors[item.pocket%4].lighter(155);edge.setAlpha(230);
                const double radius=std::clamp(item.sphere->radius*scale,2.0,100.0);
                p.setPen(QPen(edge,1.4));p.setBrush(fill);p.drawEllipse(screen.first,radius,radius);
                p.setPen(Qt::NoPen);p.setBrush(edge);p.drawEllipse(screen.first,2.0,2.0);
            }
            for(const auto& ligand:owner_->ligandCoords_){
                const auto screen=project(ligand);p.setPen(Qt::NoPen);p.setBrush(QColor(250,244,125,240));
                p.drawEllipse(screen.first,3.0,3.0);
            }
            p.setPen(QColor(230,238,248));
            p.drawText(QRect(10,8,width()-20,24),Qt::AlignLeft,
                QString("Alpha spheres: %1 | drag: rotate | wheel: zoom").arg(spheres.size()));
            return;
        }
        const bool alphaMode=owner_->algorithm_->currentIndex()==1;
        if(members.empty()){
            p.setPen(Qt::white);
            p.drawText(rect(),Qt::AlignCenter,alphaMode
                ?"Run Alpha-sphere search and select a candidate"
                :"Run Pocket Finder and select a candidate");
            return;
        }
        const double alphaSpacing=owner_->spacing_->value();
        if(alphaMode){
            QString key=QString::number(owner_->canvas_->property("pocketGeneration").toInt())
                +"|"+QString::number(alphaSpacing,'g',12);
            for(const auto id:members)key+=":"+QString::number(id);
            if(key!=alphaVoxelKey_){
                alphaVoxelPockets_.assign(owner_->alphaPockets_.size(),{});
                std::size_t totalCells=0;
                for(const auto id:members){
                    if(id>=owner_->alphaPockets_.size())continue;
                    auto& pocket=alphaVoxelPockets_[id];
                    const auto& source=owner_->alphaPockets_[id];
                    pocket.center=source.center;
                    std::set<std::array<long long,3>> cells;
                    for(const auto& sphere:source.spheres){
                        long long low[3],high[3];
                        for(int k=0;k<3;k++){
                            low[k]=static_cast<long long>(std::ceil((sphere.center[k]-sphere.radius)/alphaSpacing));
                            high[k]=static_cast<long long>(std::floor((sphere.center[k]+sphere.radius)/alphaSpacing));
                        }
                        for(long long z=low[2];z<=high[2];++z)
                            for(long long y=low[1];y<=high[1];++y)
                                for(long long x=low[0];x<=high[0];++x){
                                    const double px=x*alphaSpacing,py=y*alphaSpacing,pz=z*alphaSpacing;
                                    const double dx=px-sphere.center[0],dy=py-sphere.center[1],dz=pz-sphere.center[2];
                                    if(dx*dx+dy*dy+dz*dz>sphere.radius*sphere.radius)continue;
                                    cells.insert({x,y,z});
                                    if(cells.size()>1500000){
                                        alphaVoxelKey_.clear();
                                        p.setPen(QColor(255,220,150));
                                        p.drawText(rect(),Qt::AlignCenter,
                                            "Alpha voxel volume is too large; increase Grid spacing");
                                        return;
                                    }
                                }
                    }
                    pocket.points.reserve(cells.size());
                    for(const auto& cell:cells){
                        ligint::PocketPoint point;
                        point.pos={cell[0]*alphaSpacing,cell[1]*alphaSpacing,cell[2]*alphaSpacing};
                        point.buriedness=1.0f; // occupied alpha-sphere union, not a buriedness measurement
                        pocket.points.push_back(point);
                    }
                    pocket.pointCount=pocket.points.size();
                    pocket.volume=pocket.pointCount*alphaSpacing*alphaSpacing*alphaSpacing;
                    totalCells+=pocket.pointCount;
                    if(totalCells>3000000){
                        alphaVoxelKey_.clear();
                        p.setPen(QColor(255,220,150));
                        p.drawText(rect(),Qt::AlignCenter,
                            "Selected alpha volumes exceed the display limit; select fewer pockets or increase Grid spacing");
                        return;
                    }
                }
                alphaVoxelKey_=key;
            }
        }
        const auto& displayPockets=alphaMode?alphaVoxelPockets_:owner_->result_.pockets;
        std::vector<ligint::PocketPoint> points;
        ligint::V3 center{0,0,0};
        for(auto id:members)if(id<displayPockets.size()){
            const auto& pocket=displayPockets[id];
            points.insert(points.end(),pocket.points.begin(),pocket.points.end());
        }
        if(points.empty())return;
        // Visual-only clean-up. All operations work on a copy, never on pocket data.
        // Integer lattice indices avoid floating-point comparison during topology tests.
        const double filterSpacing=owner_->spacing_->value();
        using Cell=std::array<long long,3>;
        struct CellHash {
            std::size_t operator()(const Cell& c) const {
                std::size_t h=0;
                for(auto n:c)h^=std::hash<long long>{}(n)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);
                return h;
            }
        };
        auto lattice=[filterSpacing](const ligint::V3& p)->Cell{
            return {std::llround(p[0]/filterSpacing),
                    std::llround(p[1]/filterSpacing),
                    std::llround(p[2]/filterSpacing)};
        };
        static const Cell offsets[]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        std::unordered_map<Cell,std::size_t,CellHash> index;
        index.reserve(points.size()*2);
        for(std::size_t i=0;i<points.size();++i)index.emplace(lattice(points[i].pos),i);
        std::vector<unsigned char> keep(points.size(),1);
        const int minNeighbors=owner_->noiseNeighbors_->value();
        if(minNeighbors){
            for(std::size_t i=0;i<points.size();++i){
                const auto cell=lattice(points[i].pos);
                int neighborCount=0;
                for(const auto& offset:offsets){
                    Cell next=cell;
                    for(int k=0;k<3;k++)next[k]+=offset[k];
                    neighborCount+=index.count(next)?1:0;
                }
                if(neighborCount<minNeighbors)keep[i]=0;
            }
        }
        const int minComponent=owner_->noiseComponent_->value();
        if(minComponent>1){
            std::vector<unsigned char> visited(points.size(),0);
            for(std::size_t seed=0;seed<points.size();++seed){
                if(!keep[seed]||visited[seed])continue;
                std::queue<std::size_t> queue;
                std::vector<std::size_t> component;
                queue.push(seed);visited[seed]=1;
                while(!queue.empty()){
                    const auto i=queue.front();queue.pop();
                    component.push_back(i);
                    const auto cell=lattice(points[i].pos);
                    for(const auto& offset:offsets){
                        Cell next=cell;
                        for(int k=0;k<3;k++)next[k]+=offset[k];
                        auto found=index.find(next);
                        if(found!=index.end()&&keep[found->second]&&!visited[found->second]){
                            visited[found->second]=1;queue.push(found->second);
                        }
                    }
                }
                if(component.size()<static_cast<std::size_t>(minComponent))
                    for(auto i:component)keep[i]=0;
            }
        }
        std::unordered_set<Cell,CellHash> visible;
        visible.reserve(points.size()*2);
        std::vector<ligint::PocketPoint> cleanPoints;
        cleanPoints.reserve(points.size());
        for(std::size_t i=0;i<points.size();++i)if(keep[i]){
            visible.insert(lattice(points[i].pos));
            cleanPoints.push_back(points[i]);
        }
        const std::size_t hiddenCount=points.size()-cleanPoints.size();
        points.swap(cleanPoints);
        if(points.empty()){
            p.setPen(Qt::white);
            p.drawText(rect(),Qt::AlignCenter,"All voxels hidden by display noise filters");
            return;
        }
        for(const auto& a:points)for(int k=0;k<3;k++)center[k]+=a.pos[k]/points.size();
        double extent=1;
        for(const auto& a:points)for(int k=0;k<3;k++)extent=std::max(extent,std::fabs(a.pos[k]-center[k]));
        const double scale=0.43*std::min(width(),height())*zoom_/extent;
        struct Dot {double depth;QPointF screen;QColor color;bool protein=false;};
        struct Face {double depth;QPolygonF polygon;QColor color;};
        std::vector<Face> faces;
        std::vector<Dot> dots;dots.reserve(points.size());
        static const QColor palette[]={QColor(40,185,230),QColor(255,160,60),QColor(120,220,130),QColor(214,100,236)};
        static const QColor warm[]={QColor(248,208,77),QColor(252,107,95),QColor(151,112,244),QColor(78,196,173)};
        const double cy=std::cos(yaw_),sy=std::sin(yaw_),cp=std::cos(pitch_),sp=std::sin(pitch_);
        // Faint receptor context: local protein atoms, subsampled to keep painting responsive.
        for(std::size_t i=0;i<owner_->atoms_.size();i+=2){
            const auto& atom=owner_->atoms_[i];
            double r2=0;for(int k=0;k<3;k++)r2+=std::pow(atom.pos[k]-center[k],2);
            if(r2>400)continue;
            const double x=atom.pos[0]-center[0],y=atom.pos[1]-center[1],z=atom.pos[2]-center[2];
            const double rx=cy*x+sy*z,rz=-sy*x+cy*z;
            dots.push_back({sp*y+cp*rz,
                QPointF(width()/2.0+rx*scale,height()/2.0-(cp*y-sp*rz)*scale),
                QColor(148,159,171,75),true});
        }
        const int mode=owner_->renderMode_->currentIndex();
        const double spacing=owner_->spacing_->value();
        auto key=[spacing](const ligint::V3& v){
            return std::to_string(std::llround(v[0]/spacing))+","+
                   std::to_string(std::llround(v[1]/spacing))+","+
                   std::to_string(std::llround(v[2]/spacing));
        };
        const double dirs[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        auto project=[&](const ligint::V3& v){
            double x=v[0]-center[0],y=v[1]-center[1],z=v[2]-center[2];
            double rx=cy*x+sy*z,rz=-sy*x+cy*z;
            return std::pair<QPointF,double>{
                QPointF(width()/2.0+rx*scale,height()/2.0-(cp*y-sp*rz)*scale),
                sp*y+cp*rz};
        };
        std::unordered_set<std::string> occupied;
        if(mode!=0 && mode!=3){
            occupied.reserve(points.size()*2);
            for(const auto& voxel:points)occupied.insert(key(voxel.pos));
        }
        for(std::size_t m=0;m<members.size();m++){
            const auto id=members[m];if(id>=displayPockets.size())continue;
            for(const auto& voxel:displayPockets[id].points){
                if(!visible.count(lattice(voxel.pos)))continue;
                const double x=voxel.pos[0]-center[0],y=voxel.pos[1]-center[1],z=voxel.pos[2]-center[2];
                const double rx=cy*x+sy*z,rz=-sy*x+cy*z;
                const double ry=cp*y-sp*rz,depth=sp*y+cp*rz;
                QColor color=owner_->palette_->currentIndex()==0?palette[id%4]:warm[id%4];
                color.setAlphaF(std::clamp((owner_->opacity_->value()/100.0)*(0.3+0.7*voxel.buriedness),0.0,1.0));
                if(mode!=1 && mode!=3)
                    dots.push_back({depth,QPointF(width()/2.0+rx*scale,height()/2.0-ry*scale),color});
                if(mode==1 || mode==2){
                    for(const auto& d:dirs){
                        ligint::V3 neighbor=voxel.pos;
                        for(int k=0;k<3;k++)neighbor[k]+=d[k]*spacing;
                        if(occupied.count(key(neighbor)))continue;
                        ligint::V3 u{},v{};
                        if(d[0]!=0){u={0,spacing/2,0};v={0,0,spacing/2};}
                        else if(d[1]!=0){u={spacing/2,0,0};v={0,0,spacing/2};}
                        else{u={spacing/2,0,0};v={0,spacing/2,0};}
                        ligint::V3 midpoint=voxel.pos;
                        for(int k=0;k<3;k++)midpoint[k]+=d[k]*spacing/2;
                        QPolygonF quad;
                        double averageDepth=0;
                        for(int sign=0;sign<4;sign++){
                            const int su=(sign==0||sign==3)?-1:1;
                            const int sv=sign<2?-1:1;
                            ligint::V3 corner=midpoint;
                            for(int k=0;k<3;k++)corner[k]+=su*u[k]+sv*v[k];
                            // Consistent display-only vertex relaxation on the lattice.
                            // Shared corners see identical occupancy and therefore remain watertight.
                            if(owner_->smoothSteps_->value()>0){
                                ligint::V3 centroid{0,0,0};
                                int count=0;
                                for(int dx=-1;dx<=1;dx+=2)
                                    for(int dy=-1;dy<=1;dy+=2)
                                        for(int dz=-1;dz<=1;dz+=2){
                                            ligint::V3 sample=corner;
                                            sample[0]+=dx*spacing/2;
                                            sample[1]+=dy*spacing/2;
                                            sample[2]+=dz*spacing/2;
                                            if(!visible.count(lattice(sample)))continue;
                                            for(int k=0;k<3;k++)centroid[k]+=sample[k];
                                            ++count;
                                        }
                                if(count){
                                    const double alpha=0.17*owner_->smoothSteps_->value();
                                    for(int k=0;k<3;k++)
                                        corner[k]=(1-alpha)*corner[k]+alpha*centroid[k]/count;
                                }
                            }
                            const auto projection=project(corner);
                            quad<<projection.first;averageDepth+=projection.second/4;
                        }
                        faces.push_back({averageDepth,quad,color});
                    }
                }
            }
        }
        if(mode==3){
            QString meshKey=QString::number(owner_->canvas_->property("pocketGeneration").toInt())+"|"+QString::number(alphaMode?1:0)+"|"+QString::number(spacing,'g',12)+"|"+
                QString::number(owner_->smoothSteps_->value())+"|"+
                QString::number(owner_->noiseNeighbors_->value())+"|"+
                QString::number(owner_->noiseComponent_->value());
            for(auto id:members)meshKey+=":"+QString::number(id);
            if(meshKey!=meshKey_){
                mesh_=makeSmoothPocketMesh(points,spacing,owner_->smoothSteps_->value());
                meshKey_=meshKey;
            }
            const auto& triangles=mesh_;
            QColor shell=owner_->palette_->currentIndex()==0?QColor(64,196,224):QColor(237,160,95);
            shell.setAlphaF(owner_->opacity_->value()/100.0);
            faces.reserve(triangles.size());
            for(const auto& tri:triangles){
                QPolygonF polygon;double depth=0;
                for(const auto& vertex:tri.vertices){
                    const auto projected=project(vertex);
                    polygon<<projected.first;depth+=projected.second/3.0;
                }
                // Flat normal from triangle edges. Two-sided illumination avoids
                // artificial black faces on locally inconsistent tetrahedron winding.
                const auto& a=tri.vertices[0];
                const auto& b=tri.vertices[1];
                const auto& d=tri.vertices[2];
                const double ux=b[0]-a[0],uy=b[1]-a[1],uz=b[2]-a[2];
                const double vx=d[0]-a[0],vy=d[1]-a[1],vz=d[2]-a[2];
                double nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
                const double norm=std::sqrt(nx*nx+ny*ny+nz*nz);
                if(norm<1e-10)continue;
                nx/=norm;ny/=norm;nz/=norm;
                // Fixed camera-space key light, rotated with the viewer camera.
                const double lx=0.35,ly=0.48,lz=0.80;
                const double nrx=cy*nx+sy*nz;
                const double nrz=-sy*nx+cy*nz;
                const double nry=cp*ny-sp*nrz;
                const double ndz=sp*ny+cp*nrz;
                const double diffuse=std::fabs(nrx*lx+nry*ly+ndz*lz);
                const double strength=owner_->lightStrength_->value()/100.0;
                const double spec=std::pow(std::fabs(ndz),24.0)*owner_->shine_->value()/100.0;
                const double lum=std::clamp(0.38+strength*0.62*diffuse,0.0,1.25);
                QColor lit=QColor::fromRgb(
                    std::clamp(int(shell.red()*lum+95*spec),0,255),
                    std::clamp(int(shell.green()*lum+95*spec),0,255),
                    std::clamp(int(shell.blue()*lum+95*spec),0,255),
                    shell.alpha());
                faces.push_back({depth,polygon,lit});
            }
        }
        std::sort(dots.begin(),dots.end(),[](const Dot&a,const Dot&b){return a.depth<b.depth;});
        std::sort(faces.begin(),faces.end(),[](const Face&a,const Face&b){return a.depth<b.depth;});
        p.setPen(Qt::NoPen);
        for(const auto& f:faces){p.setBrush(f.color);p.drawPolygon(f.polygon);}
        for(const auto& atom:owner_->ligandCoords_){
            const auto screen=project(atom);
            QColor ligandColor(250,244,125);ligandColor.setAlpha(240);
            dots.push_back({screen.second,screen.first,ligandColor,false});
        }
        std::sort(dots.begin(),dots.end(),[](const Dot&a,const Dot&b){return a.depth<b.depth;});
        p.setPen(Qt::NoPen);
        for(const auto& d:dots){p.setBrush(d.color);p.drawEllipse(d.screen,d.protein?1.3:owner_->pointSize_->value()/2.0,d.protein?1.3:owner_->pointSize_->value()/2.0);}
        p.setPen(QColor(230,238,248));
        p.drawText(QRect(10,8,width()-20,24),Qt::AlignLeft,
            QString("Visible: %1 | hidden noise: %2 | drag: rotate | wheel: zoom").arg(points.size()).arg(hiddenCount));
    }
private:
    PocketWorkspace* owner_;QPoint last_;double yaw_=0.5,pitch_=0.4,zoom_=1;
    QString meshKey_;
    std::vector<SmoothTriangle> mesh_;
    QString alphaVoxelKey_;
    std::vector<ligint::PocketCandidate> alphaVoxelPockets_;
};

PocketWorkspace::PocketWorkspace(QWidget* parent):QWidget(parent){
    auto* outer=new QHBoxLayout(this);
    auto* controls=new QWidget;
    controls->setMaximumWidth(290);
    auto* column=new QVBoxLayout(controls);
    algorithm_=new QComboBox;algorithm_->setObjectName("pocket_algorithm");
    algorithm_->addItems({"Directional buriedness","Alpha spheres"});
    column->addWidget(new QLabel("Pocket search method"));
    column->addWidget(algorithm_);
    auto* form=new QFormLayout;
    spacing_=new QDoubleSpinBox;spacing_->setRange(0.5,3);spacing_->setValue(1);spacing_->setSingleStep(.25);
    buried_=new QDoubleSpinBox;buried_->setRange(.1,1);buried_->setSingleStep(.05);buried_->setValue(.5);
    ray_=new QDoubleSpinBox;ray_->setRange(2,30);ray_->setValue(10);
    surface_=new QDoubleSpinBox;surface_->setRange(1,15);surface_->setValue(5);
    probe_=new QDoubleSpinBox;probe_->setRange(.5,3);probe_->setValue(1.4);
    radius_=new QSpinBox;radius_->setRange(0,2);radius_->setValue(1);
    minPts_=new QSpinBox;minPts_->setRange(1,500);minPts_->setValue(30);
    morph_=new QComboBox;morph_->addItems({"none","opening","closing","open-close"});morph_->setCurrentText("opening");
    chain_=new QComboBox;chain_->addItem("All chains",QString());
    palette_=new QComboBox;palette_->addItems({"Cool palette","Warm palette"});
    opacity_=new QSpinBox;opacity_->setRange(10,100);opacity_->setValue(80);opacity_->setSuffix("%");
    pointSize_=new QSpinBox;pointSize_->setRange(2,12);pointSize_->setValue(5);
    // Primary selection remains visible; advanced detection parameters are collapsible.
    column->addWidget(new QLabel("Protein chain"));
    column->addWidget(chain_);
    form->addRow("Grid spacing (Å)",spacing_);form->addRow("Min buriedness",buried_);
    form->addRow("Ray length (Å)",ray_);form->addRow("Surface distance (Å)",surface_);
    form->addRow("Solvent probe (Å)",probe_);form->addRow("Morphology",morph_);
    form->addRow("Morph radius",radius_);form->addRow("Min cluster voxels",minPts_);
    renderMode_=new QComboBox;renderMode_->addItems({"Voxels","Voxel shell","Shell + Voxels","Smooth isosurface"});
    ligand_=new QComboBox;ligand_->addItem("None",QString());
    seedSpacing_=new QDoubleSpinBox;seedSpacing_->setRange(1,20);seedSpacing_->setValue(6);
    seedProminence_=new QDoubleSpinBox;seedProminence_->setRange(0,1);seedProminence_->setDecimals(2);seedProminence_->setValue(.10);
    mergeDistance_=new QDoubleSpinBox;mergeDistance_->setRange(1,30);mergeDistance_->setValue(10);
    mergeDrop_=new QDoubleSpinBox;mergeDrop_->setRange(0,1);mergeDrop_->setDecimals(2);mergeDrop_->setValue(.08);
    maxSeeds_=new QSpinBox;maxSeeds_->setRange(1,256);maxSeeds_->setValue(64);
    form->addRow("Seed separation (Å)",seedSpacing_);
    form->addRow("Seed prominence",seedProminence_);
    form->addRow("Maximum seeds",maxSeeds_);
    form->addRow("Merge max distance (Å)",mergeDistance_);
    form->addRow("Merge saddle drop",mergeDrop_);
    noiseNeighbors_=new QSpinBox;noiseNeighbors_->setRange(0,6);noiseNeighbors_->setValue(1);
    noiseComponent_=new QSpinBox;noiseComponent_->setRange(1,500);noiseComponent_->setValue(12);
    noiseNeighbors_->setToolTip("Minimum occupied face-neighbors. 0 disables this display filter.");
    noiseComponent_->setToolTip("Hide 6-connected voxel components below this size (display only).");
    form->addRow("Remove isolated voxels",noiseNeighbors_);
    form->addRow("Min visible island",noiseComponent_);
    smoothSteps_=new QSpinBox;smoothSteps_->setRange(0,4);smoothSteps_->setValue(2);
    smoothSteps_->setToolTip("Relax surface corners toward occupied voxel centers; display only. Zero retains the original voxel boundary.");
    form->addRow("Surface smoothing",smoothSteps_);
    lightStrength_=new QSpinBox;lightStrength_->setRange(0,100);lightStrength_->setValue(75);lightStrength_->setSuffix("%");
    shine_=new QSpinBox;shine_->setRange(0,100);shine_->setValue(30);shine_->setSuffix("%");
    form->addRow("Diffuse light",lightStrength_);
    form->addRow("Specular highlight",shine_);
    form->addRow("Ligand overlay",ligand_);
    form->addRow("Voxel opacity",opacity_);
    form->addRow("Voxel size",pointSize_);
    auto* alphaForm=new QFormLayout;
    alphaRender_=new QComboBox;alphaRender_->addItems({"Voxel volume","Alpha spheres"});
    alphaMinRadius_=new QDoubleSpinBox;alphaMinRadius_->setRange(0.5,10);alphaMinRadius_->setDecimals(2);alphaMinRadius_->setValue(2.0);alphaMinRadius_->setSuffix(" Å");
    alphaMaxRadius_=new QDoubleSpinBox;alphaMaxRadius_->setRange(1,15);alphaMaxRadius_->setDecimals(2);alphaMaxRadius_->setValue(5.5);alphaMaxRadius_->setSuffix(" Å");
    alphaNeighborCutoff_=new QDoubleSpinBox;alphaNeighborCutoff_->setRange(2,20);alphaNeighborCutoff_->setDecimals(1);alphaNeighborCutoff_->setValue(10);alphaNeighborCutoff_->setSuffix(" Å");
    alphaClusterDistance_=new QDoubleSpinBox;alphaClusterDistance_->setRange(0.5,10);alphaClusterDistance_->setDecimals(1);alphaClusterDistance_->setValue(3);alphaClusterDistance_->setSuffix(" Å");
    alphaMaxClusterRadius_=new QDoubleSpinBox;alphaMaxClusterRadius_->setRange(1,20);alphaMaxClusterRadius_->setDecimals(1);alphaMaxClusterRadius_->setValue(9);alphaMaxClusterRadius_->setSuffix(" Å");
    alphaMinContacts_=new QSpinBox;alphaMinContacts_->setRange(4,100);alphaMinContacts_->setValue(4);
    alphaOpacity_=new QSpinBox;alphaOpacity_->setRange(10,100);alphaOpacity_->setValue(72);alphaOpacity_->setSuffix("%");
    alphaMinSpheres_=new QSpinBox;alphaMinSpheres_->setRange(1,100);alphaMinSpheres_->setValue(3);
    alphaMacroLink_=new QDoubleSpinBox;alphaMacroLink_->setRange(0.5,15);alphaMacroLink_->setDecimals(1);alphaMacroLink_->setValue(5);alphaMacroLink_->setSuffix(" Å");
    alphaMacroDiameter_=new QDoubleSpinBox;alphaMacroDiameter_->setRange(1,40);alphaMacroDiameter_->setDecimals(1);alphaMacroDiameter_->setValue(18);alphaMacroDiameter_->setSuffix(" Å");
    alphaMacroClearance_=new QDoubleSpinBox;alphaMacroClearance_->setRange(0,5);alphaMacroClearance_->setDecimals(1);alphaMacroClearance_->setValue(1);alphaMacroClearance_->setSuffix(" Å");
    alphaMerge_=new QCheckBox("Merge connected microclusters");alphaMerge_->setChecked(true);
    alphaRank_=new QComboBox;alphaRank_->addItems({"Sphere count (default)","Sphere volume sum"});
    alphaForm->addRow("Sphere radius min",alphaMinRadius_);alphaForm->addRow("Sphere radius max",alphaMaxRadius_);
    alphaForm->addRow("Neighbor cutoff",alphaNeighborCutoff_);alphaForm->addRow("Cluster distance",alphaClusterDistance_);
    alphaForm->addRow("Max cluster radius",alphaMaxClusterRadius_);alphaForm->addRow("Min contact atoms",alphaMinContacts_);
    alphaForm->addRow("Min spheres per pocket",alphaMinSpheres_);alphaForm->addRow("Macro link distance",alphaMacroLink_);
    alphaForm->addRow("Macro max diameter",alphaMacroDiameter_);alphaForm->addRow("Macro min clearance",alphaMacroClearance_);
    alphaForm->addRow(alphaMerge_);
    alphaForm->addRow("Rank candidates by",alphaRank_);
    alphaForm->addRow("Geometry display",alphaRender_);
    alphaForm->addRow("Sphere opacity",alphaOpacity_);
    overlayMaterial_=new QComboBox;
    overlayMaterial_->addItems({"Satin translucent","Matte","Glossy","Custom"});
    overlayVisible_=new QCheckBox("Show pocket in Interactions 3D");
    overlayVisible_->setChecked(true);
    overlayOpacity_=new QSpinBox;overlayOpacity_->setRange(5,100);
    overlayOpacity_->setValue(38);overlayOpacity_->setSuffix("%");
    overlaySpecular_=new QSpinBox;overlaySpecular_->setRange(0,100);
    overlaySpecular_->setValue(12);overlaySpecular_->setSuffix("%");
    form->addRow("3D overlay opacity",overlayOpacity_);
    form->addRow("3D overlay specular",overlaySpecular_);
    auto* advancedToggle = new QToolButton;
    advancedToggle->setObjectName("pocketAdvancedToggle");
    advancedToggle->setText("Detection and display parameters");
    advancedToggle->setCheckable(true);
    advancedToggle->setChecked(false);
    advancedToggle->setArrowType(Qt::RightArrow);
    advancedToggle->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    column->addWidget(advancedToggle);
    auto* advancedArea = new QScrollArea;
    advancedArea->setObjectName("pocketAdvancedScroll");
    advancedArea->setWidgetResizable(true);
    advancedArea->setFrameShape(QFrame::StyledPanel);
    advancedArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    advancedArea->setMaximumHeight(260);
    auto* advancedContent = new QWidget;
    auto* advancedLayout = new QVBoxLayout(advancedContent);
    advancedLayout->setContentsMargins(5,5,5,5);
    auto* parameterPages=new QStackedWidget;
    auto* buriednessPage=new QWidget;buriednessPage->setLayout(form);
    auto* alphaPage=new QWidget;alphaPage->setLayout(alphaForm);
    parameterPages->addWidget(buriednessPage);parameterPages->addWidget(alphaPage);
    advancedLayout->addWidget(parameterPages);
    advancedArea->setWidget(advancedContent);
    advancedArea->setVisible(false);
    column->addWidget(advancedArea);
    connect(algorithm_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this,parameterPages,advancedToggle](int index){
        Q_UNUSED(advancedToggle);
        parameterPages->setCurrentIndex(index);
        display_->setVisible(index==0);
        run_->setText(index==0?"Find pockets":"Find Alpha pockets");
        exportAlpha_->setVisible(index==1);
        copyAlphaCommand_->setVisible(index==1);
        overlayVisible_->setVisible(index==0);
        entries_->clear();
        if(index==1) {
            for(std::size_t i=0;i<alphaPockets_.size();++i) {
                const auto& pocket=alphaPockets_[i];
                entries_->addItem(QString("#%1 | %2 spheres | %3 Å³ | center (%4, %5, %6)")
                    .arg(i+1).arg(pocket.spheres.size()).arg(pocket.volume,0,'f',0)
                    .arg(pocket.center[0],0,'f',1).arg(pocket.center[1],0,'f',1).arg(pocket.center[2],0,'f',1));
            }
        } else if(display_->currentIndex()==0) {
            for(std::size_t i=0;i<result_.pockets.size();++i)
                entries_->addItem(QString("#%1 | %2 Å³ | buried %3").arg(i+1)
                    .arg(result_.pockets[i].volume,0,'f',0).arg(result_.pockets[i].meanBuriedness,0,'f',2));
        } else {
            for(const auto& g:groups_) {
                QStringList ids;for(auto id:g.members)ids<<QString::number(id+1);
                entries_->addItem(QString("[%1] | %2 Å³ | score %3")
                    .arg(ids.join("+")).arg(g.volume,0,'f',0).arg(g.score,0,'f',3));
            }
        }
        if(entries_->count())entries_->setCurrentRow(0);
        emit pocketChoicesChanged();
        canvas_->update();publishSelection();
    });
    connect(advancedToggle,&QToolButton::toggled,this,[advancedToggle,advancedArea](bool visible){
        advancedArea->setVisible(visible);
        advancedToggle->setArrowType(visible ? Qt::DownArrow : Qt::RightArrow);
    });
    run_=new QPushButton("Find pockets");column->addWidget(run_);
    copyAlphaCommand_=new QPushButton("Copy Alpha command line");copyAlphaCommand_->setVisible(false);
    column->addWidget(copyAlphaCommand_);
    connect(copyAlphaCommand_,&QPushButton::clicked,this,[this]{
        QStringList args{path_};
        auto add=[&](const QString& key,double value){args<<key<<QString::number(value,'g',8);};
        add("--min-radius",alphaMinRadius_->value());add("--max-radius",alphaMaxRadius_->value());
        add("--neighbor-cutoff",alphaNeighborCutoff_->value());add("--cluster-distance",alphaClusterDistance_->value());
        add("--max-cluster-radius",alphaMaxClusterRadius_->value());args<<"--min-contact-atoms"<<QString::number(alphaMinContacts_->value());
        args<<"--min-spheres"<<QString::number(alphaMinSpheres_->value());
        add("--macro-link",alphaMacroLink_->value());add("--macro-diameter",alphaMacroDiameter_->value());
        add("--macro-clearance",alphaMacroClearance_->value());
        args<<"--rank-by"<<(alphaRank_->currentIndex()==1?"sphere-volume":"sphere-count");
        if(!alphaMerge_->isChecked())args<<"--no-macro-merge";
        const QString chain=chain_->currentData().toString();if(!chain.isEmpty())args<<"--chain"<<chain;
        QApplication::clipboard()->setText(joinCommand("ligint-pocket-alpha",args));
        status_->setText("Alpha command line copied.");
    });
    exportAlpha_=new QPushButton("Export Alpha results (CSV + PDB)");exportAlpha_->setEnabled(false);exportAlpha_->setVisible(false);
    exportAlpha_->setToolTip("Export candidate centers and alpha spheres, plus a merge-decision audit. This does not change the loaded fingerprint pocket.");
    column->addWidget(exportAlpha_);
    connect(exportAlpha_,&QPushButton::clicked,this,[this]{
        if(alphaPockets_.empty())return;
        const QString csvPath=QFileDialog::getSaveFileName(this,"Save Alpha candidates","alpha_pockets.csv","CSV files (*.csv)");
        if(csvPath.isEmpty())return;
        const QString stem=csvPath.endsWith(".csv",Qt::CaseInsensitive)?csvPath.left(csvPath.size()-4):csvPath;
        std::ofstream csv(stem.toStdString()+".csv");
        std::ofstream pdb(stem.toStdString()+".pdb");
        std::ofstream audit(stem.toStdString()+".audit.csv");
        std::ofstream mergeAudit(stem.toStdString()+".merge_audit.csv");
        if(!csv||!pdb||!audit||!mergeAudit){status_->setText("Could not write one or more Alpha result files.");return;}
        csv<<"rank,center_x,center_y,center_z,sphere_count,sphere_volume_sum_A3\n"<<std::fixed<<std::setprecision(3);
        audit<<"rank,sphere_count,min_radius_A,median_radius_A,max_radius_A,max_center_extent_A,bbox_diagonal_A,volume_sum_A3\n"<<std::fixed<<std::setprecision(3);
        std::size_t serial=0;
        for(std::size_t i=0;i<alphaPockets_.size();++i){
            const auto& pocket=alphaPockets_[i];
            csv<<i+1<<","<<pocket.center[0]<<","<<pocket.center[1]<<","<<pocket.center[2]<<","<<pocket.spheres.size()<<","<<pocket.volume<<"\n";
            std::vector<double> radii;double extent=0;ligint::V3 low=pocket.spheres.front().center,high=low;
            for(const auto& sphere:pocket.spheres){
                radii.push_back(sphere.radius);
                for(int k=0;k<3;++k){low[k]=std::min(low[k],sphere.center[k]);high[k]=std::max(high[k],sphere.center[k]);}
                double d2=0;for(int k=0;k<3;++k)d2+=std::pow(sphere.center[k]-pocket.center[k],2);
                extent=std::max(extent,std::sqrt(d2));
            }
            std::sort(radii.begin(),radii.end());
            double box2=0;for(int k=0;k<3;++k)box2+=std::pow(high[k]-low[k],2);
            audit<<i+1<<","<<radii.size()<<","<<radii.front()<<","<<radii[radii.size()/2]<<","<<radii.back()<<","<<extent<<","<<std::sqrt(box2)<<","<<pocket.volume<<"\n";
            for(const auto& sphere:pocket.spheres){
                if(++serial>99999)break;
                char line[160];std::snprintf(line,sizeof(line),"HETATM%5zu  C   SPH A%4zu    %8.3f%8.3f%8.3f  1.00%6.2f           C\n",serial,(i+1)%10000,sphere.center[0],sphere.center[1],sphere.center[2],sphere.radius);
                pdb<<line;
            }
        }
        pdb<<"END\n";
        mergeAudit<<"microcluster_a,microcluster_b,link_distance_A,min_clearance_A,endpoint_clearance_A,merged_center_diameter_A,accepted,reason\n"<<std::fixed<<std::setprecision(3);
        for(const auto& d:alphaAudit_.decisions)mergeAudit<<d.first<<","<<d.second<<","<<d.linkDistance<<","<<d.minClearance<<","<<d.endpointClearance<<","<<d.mergedDiameter<<","<<(d.accepted?1:0)<<","<<d.reason<<"\n";
        status_->setText(QString("Alpha results exported: %1 candidates; CSV, PDB and both audits saved.").arg(alphaPockets_.size()));
    });
    auto* exportAudit=new QPushButton("Export solvent depth audit (CSV)");
    exportAudit->setEnabled(false);
    exportAudit->setToolTip("Export per-pocket reachability, probe blockage and solvent path depth percentiles.");
    column->addWidget(exportAudit);
    connect(algorithm_,QOverload<int>::of(&QComboBox::currentIndexChanged),exportAudit,[exportAudit](int index){exportAudit->setVisible(index==0);});
    connect(exportAudit,&QPushButton::clicked,this,[this] {
        if(result_.pockets.empty())return;
        const QString path=QFileDialog::getSaveFileName(this,"Save solvent depth audit",
            "pocket_solvent_depth_audit.csv","CSV files (*.csv)");
        if(path.isEmpty())return;
        QFile file(path);
        if(!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            status_->setText("Could not write audit: "+file.errorString());
            return;
        }
        QTextStream stream(&file);
        stream<<"pocket,voxels,reachable,reachable_fraction,probe_blocked,probe_blocked_fraction,"
                 "unreachable_unblocked,unreachable_unblocked_fraction,min_A,p50_A,p90_A,max_A,"
                 "mean_A,mean_buriedness,reported_mean_A,reported_max_A,"
                 "reachable_0p8_A,reachable_fraction_0p8_A,configured_probe_radius_A,reachable_configured,reachable_fraction_configured,"
                 "reachable_2p0_A,reachable_fraction_2p0_A\n";
        auto format=[](double v){return QString::number(v,'f',4);};
        for(std::size_t i=0;i<result_.pockets.size();++i) {
            const auto& pocket=result_.pockets[i];
            const auto count=pocket.points.size();
            std::size_t blocked=0,unreachableUnblocked=0,reachable08=0,reachable20=0;
            std::vector<double> depths;
            double sum=0;
            for(const auto& pt:pocket.points) {
                reachable08 += pt.reachableProbe08;
                reachable20 += pt.reachableProbe20;
                if(pt.solventProbeBlocked)++blocked;
                if(pt.solventPathDepth>=0) {
                    depths.push_back(pt.solventPathDepth);
                    sum+=pt.solventPathDepth;
                } else if(!pt.solventProbeBlocked)++unreachableUnblocked;
            }
            std::sort(depths.begin(),depths.end());
            auto percentile=[&](double p) {
                if(depths.empty())return QString();
                const double index=p*(depths.size()-1);
                const auto lo=static_cast<std::size_t>(std::floor(index));
                const auto hi=static_cast<std::size_t>(std::ceil(index));
                return format(depths[lo]+(depths[hi]-depths[lo])*(index-lo));
            };
            const double n=count?static_cast<double>(count):1.0;
            stream<<i+1<<','<<count<<','<<depths.size()<<','<<format(depths.size()/n)<<','
                  <<blocked<<','<<format(blocked/n)<<','
                  <<unreachableUnblocked<<','<<format(unreachableUnblocked/n)<<','
                  <<percentile(0)<<','<<percentile(0.5)<<','<<percentile(0.9)<<','
                  <<percentile(1)<<','
                  <<(depths.empty()?QString():format(sum/depths.size()))<<','
                  <<format(pocket.meanBuriedness)<<','
                  <<format(pocket.meanSolventPathDepth)<<','
                  <<format(pocket.maxSolventPathDepth)<<','
                  <<reachable08<<','<<format(reachable08/n)<<','
                  <<format(probe_->value())<<','<<depths.size()<<','<<format(depths.size()/n)<<','
                  <<reachable20<<','<<format(reachable20/n)<<'\n';
        }
        status_->setText(QString("Depth audit exported: %1 pockets").arg(result_.pockets.size()));
    });
    status_=new QLabel("Open a structure to start");status_->setWordWrap(true);column->addWidget(status_);
    display_=new QComboBox;display_->addItems({"Subpockets","Candidate groups"});
    column->addWidget(display_);
    entries_=new QListWidget;entries_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    entries_->setToolTip("Ctrl-click or Shift-click to select multiple pockets or candidate groups");
    column->addWidget(entries_,1);
    auto* canvas=new PocketCanvas(this);canvas_=canvas;
    // Frequent appearance controls stay directly above the pocket viewport.
    auto* viewport = new QWidget;
    auto* viewportLayout = new QVBoxLayout(viewport);
    viewportLayout->setContentsMargins(0,0,0,0);
    viewportLayout->setSpacing(4);
    auto* quickBar = new QHBoxLayout;
    quickBar->setSpacing(6);
    quickBar->addWidget(new QLabel("Render"));
    quickBar->addWidget(renderMode_,1);
    quickBar->addWidget(new QLabel("Colors"));
    quickBar->addWidget(palette_);
    quickBar->addWidget(new QLabel("Material"));
    quickBar->addWidget(overlayMaterial_);
    quickBar->addWidget(overlayVisible_);
    viewportLayout->addLayout(quickBar);
    viewportLayout->addWidget(canvas,1);
    outer->addWidget(controls);
    outer->addWidget(viewport,1);
    connect(display_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this]{
        entries_->clear();
        if(algorithm_->currentIndex()==1) {
            for(std::size_t i=0;i<alphaPockets_.size();++i) {
                const auto& pocket=alphaPockets_[i];
                entries_->addItem(QString("#%1 | %2 spheres | %3 Å³ | center (%4, %5, %6)")
                    .arg(i+1).arg(pocket.spheres.size()).arg(pocket.volume,0,'f',0)
                    .arg(pocket.center[0],0,'f',1).arg(pocket.center[1],0,'f',1).arg(pocket.center[2],0,'f',1));
            }
        } else if(display_->currentIndex()==0){
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
        emit pocketChoicesChanged();
        canvas_->update();
        publishSelection();
    });
    connect(entries_,&QListWidget::itemSelectionChanged,canvas_,QOverload<>::of(&QWidget::update));
    connect(entries_,&QListWidget::itemSelectionChanged,this,&PocketWorkspace::publishSelection);
    connect(palette_,QOverload<int>::of(&QComboBox::currentIndexChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(renderMode_,QOverload<int>::of(&QComboBox::currentIndexChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(alphaRender_,QOverload<int>::of(&QComboBox::currentIndexChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(ligand_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this]{
        ligandCoords_.clear();
        const QString token=ligand_->currentData().toString();
        if(!token.isEmpty()){
            std::ifstream file(path_.toStdString());
            std::string line;
            while(std::getline(file,line)){
                if(line.size()<54||line.compare(0,6,"HETATM")!=0)continue;
                const QString code=QString::fromStdString(line.substr(17,3)+"|"+line.substr(21,1)+"|"+line.substr(22,4));
                if(code!=token)continue;
                if(line.size()>=78){
                    const QString el=QString::fromStdString(line.substr(76,2)).trimmed().toUpper();
                    if(el=="H"||el=="D")continue;
                }
                try{ligandCoords_.push_back({std::stod(line.substr(30,8)),
                    std::stod(line.substr(38,8)),std::stod(line.substr(46,8))});}
                catch(const std::exception&){}
            }
        }
        canvas_->update();
    });
    connect(opacity_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(alphaOpacity_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(pointSize_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(noiseNeighbors_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(noiseComponent_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(smoothSteps_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(smoothSteps_,QOverload<int>::of(&QSpinBox::valueChanged),this,[this]{publishSelection();});
    connect(lightStrength_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(shine_,QOverload<int>::of(&QSpinBox::valueChanged),canvas_,QOverload<>::of(&QWidget::update));
    connect(overlayMaterial_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this](int preset){
        if(preset==3)return;
        const QSignalBlocker blockOpacity(overlayOpacity_);
        const QSignalBlocker blockSpecular(overlaySpecular_);
        if(preset==0){overlayOpacity_->setValue(38);overlaySpecular_->setValue(12);}
        if(preset==1){overlayOpacity_->setValue(65);overlaySpecular_->setValue(0);}
        if(preset==2){overlayOpacity_->setValue(55);overlaySpecular_->setValue(55);}
        emit pocketOpacityPercentChanged(overlayOpacity_->value());
        emit pocketAppearanceChanged(overlayVisible_->isChecked(),
                                     overlayOpacity_->value()/100.0,
                                     overlaySpecular_->value()/100.0);
    });
    auto appearanceChanged=[this]{
        emit pocketAppearanceChanged(overlayVisible_->isChecked(),
                                     overlayOpacity_->value()/100.0,
                                     overlaySpecular_->value()/100.0);
    };
    connect(overlayVisible_,&QCheckBox::toggled,this,appearanceChanged);
    connect(overlayVisible_,&QCheckBox::toggled,this,&PocketWorkspace::pocketVisibilityChanged);
    connect(overlayMaterial_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,&PocketWorkspace::pocketMaterialChanged);
    connect(overlayOpacity_,QOverload<int>::of(&QSpinBox::valueChanged),this,&PocketWorkspace::pocketOpacityPercentChanged);
    connect(overlayOpacity_,QOverload<int>::of(&QSpinBox::valueChanged),this,[this,appearanceChanged](int){
        if(overlayMaterial_->currentIndex()!=3)overlayMaterial_->setCurrentIndex(3);
        appearanceChanged();
    });
    connect(overlaySpecular_,QOverload<int>::of(&QSpinBox::valueChanged),this,[this,appearanceChanged](int){
        if(overlayMaterial_->currentIndex()!=3)overlayMaterial_->setCurrentIndex(3);
        appearanceChanged();
    });
    connect(run_,&QPushButton::clicked,this,[this]{
        if(path_.isEmpty())return;
        const QString path=path_;
        const QString chain=chain_->currentData().toString();
        if(algorithm_->currentIndex()==1) {
            ligint::AlphaPocketOptions options;
            options.minRadius=alphaMinRadius_->value();options.maxRadius=alphaMaxRadius_->value();
            options.neighborCutoff=alphaNeighborCutoff_->value();options.clusterDistance=alphaClusterDistance_->value();
            options.maxClusterRadius=alphaMaxClusterRadius_->value();
            options.minContactAtoms=static_cast<std::size_t>(alphaMinContacts_->value());
            options.minSpheres=static_cast<std::size_t>(alphaMinSpheres_->value());
            options.macroLinkDistance=alphaMacroLink_->value();options.macroMaxDiameter=alphaMacroDiameter_->value();
            options.macroMinClearance=alphaMacroClearance_->value();options.mergeMicroclusters=alphaMerge_->isChecked();
            options.ranking=alphaRank_->currentIndex()==1?ligint::AlphaPocketRanking::SphereVolume:ligint::AlphaPocketRanking::SphereCount;
            run_->setEnabled(false);status_->setText("Computing Alpha-sphere pockets...");
            watcher_.setFuture(QtConcurrent::run([path,chain,options]{
                PocketGuiOutcome out;out.alphaMode=true;
                try{
                    out.atoms=ligint::readPocketAtoms(path.toStdString());
                    auto selected=out.atoms;
                    if(!chain.isEmpty())selected.erase(std::remove_if(selected.begin(),selected.end(),[&](const auto& atom){
                        return QString::fromStdString(atom.chain)!=chain;
                    }),selected.end());
                    out.alphaPockets=ligint::findAlphaPockets(selected,options,&out.alphaAudit);
                }catch(const std::exception& e){out.error=QString::fromUtf8(e.what());}
                return out;
            }));
            return;
        }
        ligint::PocketOptions o;
        o.spacing=spacing_->value();o.minBuriedness=buried_->value();o.rayLength=ray_->value();
        o.maxSurfaceDistance=surface_->value();o.solventProbeRadius=probe_->value();
        o.morphology=morph_->currentText().toStdString();o.morphologyRadius=radius_->value();
        o.minClusterPoints=static_cast<std::size_t>(minPts_->value());
        o.seedSpacing=seedSpacing_->value();o.seedProminence=seedProminence_->value();
        o.maxSeeds=static_cast<std::size_t>(maxSeeds_->value());
        o.mergeMaxSeedDistance=mergeDistance_->value();
        o.mergeSaddleDrop=mergeDrop_->value();
        if(!chain.isEmpty())o.chains={chain.toStdString()};
        run_->setEnabled(false);status_->setText("Computing pockets...");
        watcher_.setFuture(QtConcurrent::run([path,o]{
            PocketGuiOutcome out;
            try{
                auto atoms=ligint::readPocketAtoms(path.toStdString());
                out.atoms=atoms;
                out.result=ligint::findPockets(atoms,o);
                std::stable_sort(out.result.pockets.begin(),out.result.pockets.end(),
                    [](const auto&a,const auto&b){return a.heuristicScore>b.heuristicScore;});
                const auto edges=ligint::analyzePocketInterfaces(out.result.pockets,atoms,o.spacing);
                out.groups=ligint::rankPocketGroups(out.result.pockets,edges);
            }catch(const std::exception& e){out.error=QString::fromUtf8(e.what());}
            return out;
        }));
    });
    connect(&watcher_,&QFutureWatcher<PocketGuiOutcome>::finished,this,[this,exportAudit]{
        auto out=watcher_.result();run_->setEnabled(true);
        if(!out.error.isEmpty()){status_->setText("Error: "+out.error);return;}
        atoms_=std::move(out.atoms);
        if(out.alphaMode) {
            alphaPockets_=std::move(out.alphaPockets);alphaAudit_=std::move(out.alphaAudit);
            exportAlpha_->setEnabled(!alphaPockets_.empty());
            status_->setText(QString("%1 Alpha pockets; %2 spheres; %3 atoms")
                .arg(alphaPockets_.size()).arg([this]{std::size_t n=0;for(const auto& p:alphaPockets_)n+=p.spheres.size();return n;}())
                .arg(atoms_.size()));
        } else {
            result_=std::move(out.result);groups_=std::move(out.groups);
            exportAudit->setEnabled(!result_.pockets.empty());
            status_->setText(QString("%1 pockets, %2 groups; %3 atoms")
                .arg(result_.pockets.size()).arg(groups_.size()).arg(result_.selectedAtoms));
        }
        canvas_->setProperty("pocketGeneration",canvas_->property("pocketGeneration").toInt()+1);
        if(algorithm_->currentIndex()==1) {
            entries_->clear();
            for(std::size_t i=0;i<alphaPockets_.size();++i) {
                const auto& pocket=alphaPockets_[i];
                entries_->addItem(QString("#%1 | %2 spheres | %3 Å³ | center (%4, %5, %6)")
                    .arg(i+1).arg(pocket.spheres.size()).arg(pocket.volume,0,'f',0)
                    .arg(pocket.center[0],0,'f',1).arg(pocket.center[1],0,'f',1).arg(pocket.center[2],0,'f',1));
            }
        } else {
            const int mode=display_->currentIndex();
            display_->setCurrentIndex(mode==0?1:0);display_->setCurrentIndex(mode);
        }
        if(entries_->count())entries_->setCurrentRow(0);
        emit pocketChoicesChanged();canvas_->update();publishSelection();
    });
}
PocketWorkspace::~PocketWorkspace(){if(watcher_.isRunning())watcher_.waitForFinished();}
void PocketWorkspace::setStructurePath(const QString& path){
    path_=path;run_->setEnabled(!path.isEmpty());
    status_->setText(path.isEmpty()?"Open a structure to start":"Ready: "+path);
    result_={};groups_.clear();alphaPockets_.clear();alphaAudit_={};atoms_.clear();entries_->clear();exportAlpha_->setEnabled(false);
    emit pocketChoicesChanged();
    canvas_->setProperty("pocketGeneration",canvas_->property("pocketGeneration").toInt()+1);
    chain_->clear();chain_->addItem("All chains",QString());
    ligand_->clear();ligand_->addItem("None",QString());ligandCoords_.clear();
    if(!path.isEmpty()){
        try{
            const auto atoms=ligint::readPocketAtoms(path.toStdString());
            std::set<std::string> chains;
            for(const auto& atom:atoms)if(!atom.chain.empty())chains.insert(atom.chain);
            for(const auto& chain:chains)chain_->addItem(QString::fromStdString(chain),QString::fromStdString(chain));
            std::ifstream input(path.toStdString());
            std::set<std::string> ligandIDs;
            std::string line;
            while(std::getline(input,line)){
                if(line.size()<54||line.compare(0,6,"HETATM")!=0)continue;
                const std::string token=line.substr(17,3)+"|"+line.substr(21,1)+"|"+line.substr(22,4);
                if(token.substr(0,3)=="HOH")continue;
                ligandIDs.insert(token);
            }
            for(const auto& token:ligandIDs)
                ligand_->addItem(QString::fromStdString(token),QString::fromStdString(token));
        }catch(const std::exception& e){status_->setText(QString::fromUtf8(e.what()));}
    }
    canvas_->update();
    publishSelection();
}
std::vector<std::size_t> PocketWorkspace::visibleMembers() const{
    std::set<std::size_t> ids;
    for(const auto* item:entries_->selectedItems()){
        const int row=entries_->row(item);
        if(row<0)continue;
        if(algorithm_->currentIndex()==1){
            if(static_cast<std::size_t>(row)<alphaPockets_.size())ids.insert(static_cast<std::size_t>(row));
        }else if(display_->currentIndex()==0){
            if(static_cast<std::size_t>(row)<result_.pockets.size())
                ids.insert(static_cast<std::size_t>(row));
        }else if(static_cast<std::size_t>(row)<groups_.size()){
            for(auto id:groups_[row].members)ids.insert(id);
        }
    }
    return {ids.begin(),ids.end()};
}
bool PocketWorkspace::pocketOverlayVisible() const { return overlayVisible_->isChecked(); }
void PocketWorkspace::setPocketOverlayVisible(bool visible) { overlayVisible_->setChecked(visible); }
int PocketWorkspace::pocketMaterialIndex() const { return overlayMaterial_->currentIndex(); }
void PocketWorkspace::setPocketMaterialIndex(int index) {
    if (index >= 0 && index < overlayMaterial_->count()) overlayMaterial_->setCurrentIndex(index);
}
QList<int> PocketWorkspace::selectedPocketRows() const {
    QList<int> rows;
    for (auto* item : entries_->selectedItems()) rows.append(entries_->row(item));
    std::sort(rows.begin(),rows.end());
    return rows;
}
void PocketWorkspace::setSelectedPocketRows(const QList<int>& rows) {
    const QSignalBlocker blocker(entries_);
    entries_->clearSelection();
    for (int row : rows)
        if (row >= 0 && row < entries_->count())
            entries_->item(row)->setSelected(true);
    canvas_->update();
    publishSelection();
}
int PocketWorkspace::pocketOpacityPercent() const { return overlayOpacity_->value(); }
void PocketWorkspace::setPocketOpacityPercent(int percent) { overlayOpacity_->setValue(percent); }
std::vector<ligint::PocketPoint> PocketWorkspace::pocketColorSample(int row) const {
    if (row < 0 || algorithm_->currentIndex()==1) return {};
    if (display_->currentIndex() == 0) {
        if (static_cast<std::size_t>(row) < result_.pockets.size())
            { const auto& pts = result_.pockets[static_cast<std::size_t>(row)].points;
              return pts.empty() ? std::vector<ligint::PocketPoint>{} : std::vector<ligint::PocketPoint>{pts.front()}; }
    } else if (static_cast<std::size_t>(row) < groups_.size()) {
        const auto& members = groups_[static_cast<std::size_t>(row)].members;
        if (!members.empty() && members.front() < result_.pockets.size())
            { const auto& pts = result_.pockets[members.front()].points;
              return pts.empty() ? std::vector<ligint::PocketPoint>{} : std::vector<ligint::PocketPoint>{pts.front()}; }
    }
    return {};
}
QStringList PocketWorkspace::pocketChoices() const {
    QStringList choices;
    for (int i = 0; i < entries_->count(); ++i) choices << entries_->item(i)->text();
    return choices;
}
int PocketWorkspace::rowForVisiblePocket(int visibleIndex) const {
    const auto members = visibleMembers();
    if (visibleIndex < 0 || static_cast<std::size_t>(visibleIndex) >= members.size()) return -1;
    const auto id = members[static_cast<std::size_t>(visibleIndex)];
    if (algorithm_->currentIndex()==1)
        return id < static_cast<std::size_t>(entries_->count()) ? static_cast<int>(id) : -1;
    if (display_->currentIndex() == 0)
        return id < static_cast<std::size_t>(entries_->count()) ? static_cast<int>(id) : -1;
    for (int row = 0; row < entries_->count() && static_cast<std::size_t>(row) < groups_.size(); ++row) {
        const auto& group = groups_[static_cast<std::size_t>(row)];
        if (std::find(group.members.begin(), group.members.end(), id) != group.members.end())
            return row;
    }
    return -1;
}
int PocketWorkspace::selectedPocketRow() const {
    return entries_->currentRow();
}
void PocketWorkspace::selectPocketRow(int row) {
    if (row < 0 || row >= entries_->count()) return;
    entries_->setCurrentRow(row);
    canvas_->update();
    publishSelection();
}
void PocketWorkspace::publishSelection(){
    std::vector<ligint::PocketPoint> points;
    std::vector<std::vector<ligint::PocketPoint>> groups;
    if(algorithm_->currentIndex()!=1)for(const auto id:visibleMembers())if(id<result_.pockets.size()){
        const auto& source=result_.pockets[id].points;
        points.insert(points.end(),source.begin(),source.end());
        groups.push_back(source);
    }
    emit pocketOverlayChanged(points,spacing_->value(),smoothSteps_->value(),groups);
    emit pocketSelectionChanged(entries_->currentRow());
    emit pocketSelectedRowsChanged(selectedPocketRows());
}
} // namespace gui
