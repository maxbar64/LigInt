#include "Molecular3DWidget.h"

#include <QCheckBox>
#include <QComboBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QUrl>
#include <QVBoxLayout>
#include <QWebEnginePage>
#include <QWebEngineSettings>
#include <QWebEngineView>

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

namespace gui {
namespace {

QString viewerHtml() {
    return QStringLiteral(R"HTML(
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
html,body,#viewport{width:100%;height:100%;margin:0;overflow:hidden;background:#101720;color:#d9e2ec;font:14px sans-serif}
#status{position:absolute;z-index:2;left:14px;top:12px;padding:7px 10px;background:#172331cc;border-radius:6px;max-width:80%}
</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js"
        onerror="document.getElementById('status').textContent='3Dmol.js is unreachable: check the Internet connection.'"></script>
</head>
<body>
<div id="viewport"></div><div id="status">Loading the WebGL renderer…</div>
<script>
(function(){
  const status = document.getElementById('status');
  let viewer = null, model = null, ligandSelection = null, currentData = null, ligandFound = false;
  let styleIndex = 0, showContacts = true, showSurface = false;

  function statusText(text){ status.textContent = text; }
  function contactColor(type) {
    const t = String(type || '').toLowerCase();
    if (t.includes('clash') || t.includes('repulsion')) return '#ff5c5c';
    if (t.includes('hydrogen') || t.includes('metal')) return '#32d5ff';
    if (t.includes('salt') || t.includes('cation')) return '#d47bff';
    if (t.includes('pi-stacking')) return '#ffc857';
    return '#62d98b';
  }
  // Redraws styles, contacts and surface. The camera is only moved when `zoom` is true,
  // so toggling an option does not throw away the user's current zoom and pan.
  function draw(zoom) {
    if (!viewer || !model || !currentData || !ligandFound) return;
    viewer.removeAllSurfaces();
    viewer.removeAllShapes();
    const proteinStyle = styleIndex === 1
      ? {line:{color:'#b7c4d4'}}
      : styleIndex === 2
        ? {cartoon:{color:'#c4d0df',ribbon:true,thickness:0.5}}
        : {cartoon:{color:'spectrum'}};
    viewer.setStyle({hetflag:false}, proteinStyle);
    viewer.setStyle(ligandSelection, {
      stick:{colorscheme:'Jmol',radius:0.2},
      sphere:{scale:0.23,colorscheme:'Jmol'}
    });
    if (showContacts) {
      (currentData.hits || []).forEach(function(hit) {
        const a = hit.ligand, b = hit.protein;
        viewer.addLine({
          start:{x:a[0],y:a[1],z:a[2]},
          end:{x:b[0],y:b[1],z:b[2]},
          color:contactColor(hit.type), linewidth:2, dashed:true
        });
      });
    }
    if (showSurface && window.$3Dmol) {
      viewer.addSurface($3Dmol.SurfaceType.SAS, {
        opacity:0.62, color:'#58a6d8'
      }, {within:{distance:5,sel:ligandSelection},hetflag:false}, {hetflag:false}, ligandSelection);
    }
    if (zoom) viewer.zoomTo(ligandSelection);
    viewer.render();
  }

  window.setLigIntStyle = function(index){ styleIndex = Number(index) || 0; draw(false); };
  window.setLigIntContacts = function(enabled){ showContacts = !!enabled; draw(false); };
  window.setLigIntSurface = function(enabled){ showSurface = !!enabled; draw(false); };
  window.fitLigIntView = function(){ if(viewer && ligandFound) { viewer.zoomTo(ligandSelection); viewer.render(); } };
  window.fitLigIntProtein = function(){ if(viewer && model) { viewer.zoomTo(); viewer.render(); } };
  window.clearLigInt = function(){
    model = null; currentData = null; ligandFound = false;
    if (viewer) { viewer.clear(); viewer.render(); }
    statusText('No structure loaded.');
  };
  window.renderLigInt = function(data) {
    if (!window.$3Dmol) {
      statusText('3Dmol.js could not be loaded. Check the Internet connection.');
      return;
    }
    currentData = data;
    const opts = data.options || {};
    styleIndex = Number(opts.style) || 0;
    showContacts = opts.contacts !== false;
    showSurface = !!opts.surface;
    // One viewer for the lifetime of the page: creating a new one on every call would leak
    // a canvas and a WebGL context each time a different ligand is analyzed.
    if (!viewer) viewer = $3Dmol.createViewer(document.getElementById('viewport'), {backgroundColor:'#101720', antialias:true});
    else viewer.clear();
    model = viewer.addModel(data.structure, data.format);
    const lig = data.ligand || {};
    ligandSelection = {resn:lig.name, resi:Number(lig.resi)};
    if (lig.chain) ligandSelection.chain = lig.chain;
    if (lig.icode && lig.icode.trim()) ligandSelection.icode = lig.icode;
    ligandFound = model.selectedAtoms(ligandSelection).length > 0;
    if (!ligandFound) {
      const relaxed = {resn:lig.name};
      if (lig.chain) relaxed.chain = lig.chain;
      ligandSelection = relaxed;
      ligandFound = model.selectedAtoms(ligandSelection).length > 0;
    }
    if (!ligandFound) {
      statusText('Protein loaded, but the selected ligand could not be located.');
      viewer.setStyle({hetflag:false},{cartoon:{color:'spectrum'}});
      viewer.zoomTo();
      viewer.render();
      return;
    }
    statusText('Whole protein · ligand ball-and-stick · analysis contacts');
    draw(true);
    setTimeout(function(){ if(viewer) viewer.resize(); }, 100);
  };
  window.addEventListener('resize', function(){ if(viewer) viewer.resize(); });
})();
</script>
</body>
</html>)HTML");
}

}  // namespace

Molecular3DWidget::Molecular3DWidget(QWidget* parent) : QWidget(parent) {
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(6, 6, 6, 6);
    layout->setSpacing(6);

    auto* controls = new QHBoxLayout;
    controls->addWidget(new QLabel(tr("Protein style:"), this));
    styleCombo_ = new QComboBox(this);
    styleCombo_->setObjectName("advanced3d_protein_style");
    styleCombo_->addItem(tr("Cartoon"));
    styleCombo_->addItem(tr("Line"));
    styleCombo_->addItem(tr("Soft ribbon"));
    controls->addWidget(styleCombo_);

    contactsCheck_ = new QCheckBox(tr("Interaction contacts"), this);
    contactsCheck_->setObjectName("advanced3d_contacts");
    contactsCheck_->setChecked(true);
    controls->addWidget(contactsCheck_);

    surfaceCheck_ = new QCheckBox(tr("Pocket surface"), this);
    surfaceCheck_->setObjectName("advanced3d_surface");
    controls->addWidget(surfaceCheck_);

    auto* fit = new QPushButton(tr("Focus ligand"), this);
    fit->setObjectName("advanced3d_focus");
    controls->addWidget(fit);
    auto* fitProtein = new QPushButton(tr("Fit protein"), this);
    fitProtein->setObjectName("advanced3d_fit_protein");
    controls->addWidget(fitProtein);
    layout->addLayout(controls);

    webView_ = new QWebEngineView(this);
    webView_->setObjectName("advanced3d_view");
    webView_->settings()->setAttribute(QWebEngineSettings::LocalContentCanAccessRemoteUrls, true);
    layout->addWidget(webView_, 1);

    connect(webView_, &QWebEngineView::loadFinished, this, [this](bool ok) {
        pageReady_ = ok;
        if (ok) sendPayload();
    });
    connect(styleCombo_, qOverload<int>(&QComboBox::currentIndexChanged), this,
            [this](int index) { runViewerCommand(QStringLiteral("setLigIntStyle(%1)").arg(index)); });
    connect(contactsCheck_, &QCheckBox::toggled, this,
            [this](bool enabled) { runViewerCommand(QStringLiteral("setLigIntContacts(%1)").arg(enabled ? "true" : "false")); });
    connect(surfaceCheck_, &QCheckBox::toggled, this,
            [this](bool enabled) { runViewerCommand(QStringLiteral("setLigIntSurface(%1)").arg(enabled ? "true" : "false")); });
    connect(fit, &QPushButton::clicked, this, [this] { runViewerCommand(QStringLiteral("fitLigIntView()")); });
    connect(fitProtein, &QPushButton::clicked, this, [this] { runViewerCommand(QStringLiteral("fitLigIntProtein()")); });

    webView_->setHtml(viewerHtml(), QUrl(QStringLiteral("https://ligint.local/")));
}

void Molecular3DWidget::setStructure(const QString& data, const QString& format,
                                     const QString& ligandName, const QString& chain,
                                     int residueNumber, char insertionCode,
                                     const std::vector<ligint::HitRecord>& hits) {
    QJsonObject ligand;
    ligand.insert(QStringLiteral("name"), ligandName);
    ligand.insert(QStringLiteral("chain"), chain);
    ligand.insert(QStringLiteral("resi"), residueNumber);
    // A blank or NUL insertion code means "none": sending it would make the JS selection match nothing.
    ligand.insert(QStringLiteral("icode"),
                  (insertionCode == ' ' || insertionCode == '\0') ? QString() : QString(QChar::fromLatin1(insertionCode)));
    payload_.insert(QStringLiteral("structure"), data);
    payload_.insert(QStringLiteral("format"), format);
    payload_.insert(QStringLiteral("ligand"), ligand);

    QJsonArray hitArray;
    for (const auto& hit : hits) {
        QJsonObject item;
        item.insert(QStringLiteral("type"), QString::fromStdString(hit.type));
        QJsonArray ligandPos, proteinPos;
        for (double v : hit.ligandPos) ligandPos.append(v);
        for (double v : hit.proteinPos) proteinPos.append(v);
        item.insert(QStringLiteral("ligand"), ligandPos);
        item.insert(QStringLiteral("protein"), proteinPos);
        hitArray.append(item);
    }
    payload_.insert(QStringLiteral("hits"), hitArray);
    sendPayload();
}

void Molecular3DWidget::clear() {
    payload_ = QJsonObject();
    runViewerCommand(QStringLiteral("clearLigInt()"));
}

void Molecular3DWidget::sendPayload() {
    if (!pageReady_ || payload_.isEmpty()) return;
    // The current control state travels with the payload, so the page draws the scene once
    // instead of once per option.
    QJsonObject options;
    options.insert(QStringLiteral("style"), styleCombo_->currentIndex());
    options.insert(QStringLiteral("contacts"), contactsCheck_->isChecked());
    options.insert(QStringLiteral("surface"), surfaceCheck_->isChecked());
    QJsonObject payload = payload_;
    payload.insert(QStringLiteral("options"), options);
    const QString json = QString::fromUtf8(QJsonDocument(payload).toJson(QJsonDocument::Compact));
    webView_->page()->runJavaScript(QStringLiteral("window.renderLigInt(%1);").arg(json));
}

void Molecular3DWidget::runViewerCommand(const QString& command) {
    if (pageReady_) webView_->page()->runJavaScript(QStringLiteral("window.%1;").arg(command));
}

}  // namespace gui
