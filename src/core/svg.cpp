// SVG figure of an interaction fingerprint (residues grouped by pocket region, one glyph per interaction class).
#include <algorithm>
#include <cstdio>
#include <sstream>

#include "ligint/fingerprint.hpp"

namespace ligint {
namespace {

std::string esc(const std::string& s) {
    std::string r;
    for (char c : s) {
        switch (c) {
            case '&': r += "&amp;"; break;
            case '<': r += "&lt;"; break;
            case '>': r += "&gt;"; break;
            case '"': r += "&quot;"; break;
            default: r += c;
        }
    }
    return r;
}

std::string f1(double d) {
    char b[32];
    std::snprintf(b, sizeof b, "%.1f", d);
    return b;
}

// color and glyph of each class
const char* classColor(int c) {
    static const char* col[FpClassCount] = {"#868e96", "#2f9e44", "#94d82d", "#e03131", "#1c7ed6", "#1971c2", "#c92a2a",
                                            "#f08c00", "#0ca678", "#7048e8", "#4dabf7"};
    return col[c];
}

std::string glyph(int c, double cx, double cy, double r) {
    std::ostringstream s;
    const std::string col = classColor(c);
    auto pt = [&](double x, double y) { return f1(x) + "," + f1(y); };
    switch (c) {
        case FpHydrophobic:          // circle
        case FpFaceToFace:
        case FpEdgeToFace:
            s << "<circle cx=\"" << f1(cx) << "\" cy=\"" << f1(cy) << "\" r=\"" << f1(r) << "\" fill=\"" << col << "\"/>";
            break;
        case FpHBondDonor:           // triangle pointing up
            s << "<polygon points=\"" << pt(cx, cy - r) << " " << pt(cx + r, cy + r) << " " << pt(cx - r, cy + r) << "\" fill=\"" << col << "\"/>";
            break;
        case FpHBondAcceptor:        // triangle pointing down
            s << "<polygon points=\"" << pt(cx, cy + r) << " " << pt(cx + r, cy - r) << " " << pt(cx - r, cy - r) << "\" fill=\"" << col << "\"/>";
            break;
        case FpIonicPositive:        // circle with a plus
            s << "<circle cx=\"" << f1(cx) << "\" cy=\"" << f1(cy) << "\" r=\"" << f1(r) << "\" fill=\"" << col << "\"/>"
              << "<path d=\"M" << pt(cx - r * 0.55, cy) << " L" << pt(cx + r * 0.55, cy) << " M" << pt(cx, cy - r * 0.55) << " L" << pt(cx, cy + r * 0.55)
              << "\" stroke=\"#fff\" stroke-width=\"0.9\"/>";
            break;
        case FpIonicNegative:        // circle with a minus
            s << "<circle cx=\"" << f1(cx) << "\" cy=\"" << f1(cy) << "\" r=\"" << f1(r) << "\" fill=\"" << col << "\"/>"
              << "<path d=\"M" << pt(cx - r * 0.55, cy) << " L" << pt(cx + r * 0.55, cy) << "\" stroke=\"#fff\" stroke-width=\"0.9\"/>";
            break;
        case FpCationPi:             // diamond
            s << "<polygon points=\"" << pt(cx, cy - r) << " " << pt(cx + r, cy) << " " << pt(cx, cy + r) << " " << pt(cx - r, cy) << "\" fill=\"" << col << "\"/>";
            break;
        case FpHalogen:              // square
        case FpMetal:
            s << "<rect x=\"" << f1(cx - r * 0.85) << "\" y=\"" << f1(cy - r * 0.85) << "\" width=\"" << f1(r * 1.7) << "\" height=\"" << f1(r * 1.7)
              << "\" fill=\"" << col << "\"/>";
            break;
        default:                     // water bridge: hollow circle
            s << "<circle cx=\"" << f1(cx) << "\" cy=\"" << f1(cy) << "\" r=\"" << f1(r * 0.85) << "\" fill=\"none\" stroke=\"" << col << "\" stroke-width=\"1.3\"/>";
    }
    return s.str();
}

// color of the residue letter: the most informative class present
const char* letterColor(const FingerprintEntry& e) {
    static const int prio[] = {FpIonicPositive, FpIonicNegative, FpHBondDonor, FpHBondAcceptor, FpFaceToFace, FpEdgeToFace,
                               FpCationPi, FpHalogen, FpMetal, FpWaterBridge, FpHydrophobic};
    for (int c : prio)
        if (e.counts[c] > 0) return classColor(c);
    return "#adb5bd";
}

struct Chunk {
    std::string region;
    std::vector<const FingerprintEntry*> cells;
};

}  // namespace

std::string fingerprintSvg(const Fingerprint& fp, const SvgOptions& opt) {
    const int cellsPerRow = std::max(4, opt.cellsPerRow);
    const double cw = 32, margin = 12, headerH = 15, labelH = 11, letterH = 14, numberH = 11, glyphH = 27, gap = 14;
    const double cellH = labelH + letterH + numberH + glyphH;

    // ---- regions, in order of appearance
    auto regionOf = [&](const FingerprintEntry& e) -> std::string {
        if (!e.inPocket) return "other contacts";
        return e.region.empty() ? "contacts" : e.region;
    };
    std::vector<std::vector<Chunk>> rows;
    {
        int remaining = 0;
        size_t i = 0;
        while (i < fp.entries.size()) {
            const std::string reg = regionOf(fp.entries[i]);
            size_t j = i;
            while (j < fp.entries.size() && regionOf(fp.entries[j]) == reg) j++;
            size_t idx = i;
            // keep a region whole when it fits in a row: start a new row instead of splitting it
            if (j - i <= static_cast<size_t>(cellsPerRow) && static_cast<size_t>(remaining) < j - i) remaining = 0;
            while (idx < j) {
                if (remaining == 0) { rows.emplace_back(); remaining = cellsPerRow; }
                const size_t n = std::min<size_t>(remaining, j - idx);
                Chunk c;
                c.region = reg;
                for (size_t k = 0; k < n; k++) c.cells.push_back(&fp.entries[idx + k]);
                rows.back().push_back(c);
                remaining -= static_cast<int>(n);
                idx += n;
            }
            i = j;
        }
    }
    size_t maxCells = 0;
    for (const auto& r : rows) {
        size_t n = 0;
        for (const auto& c : r) n += c.cells.size();
        maxCells = std::max(maxCells, n);
    }

    // region colors: from the pocket definition when given, otherwise a palette in order of first appearance
    static const char* palette[] = {"#f3e04a", "#7bc96f", "#e5857f", "#f4a261", "#6fa8dc", "#c49bd6", "#f5b041", "#76d7c4",
                                    "#bfc9d1", "#ee9ca7", "#a3c4f3", "#d4e09b"};
    std::vector<std::string> regionOrder;
    for (const auto& r : rows)
        for (const auto& c : r)
            if (std::find(regionOrder.begin(), regionOrder.end(), c.region) == regionOrder.end()) regionOrder.push_back(c.region);
    auto colorOfRegion = [&](const std::string& reg) {
        for (const auto& rc : fp.regionColors)
            if (rc.first == reg) return rc.second;
        size_t k = std::find(regionOrder.begin(), regionOrder.end(), reg) - regionOrder.begin();
        return std::string(palette[k % (sizeof palette / sizeof palette[0])]);
    };

    // ---- legend layout
    const double width = 2 * margin + std::max<double>(maxCells * cw, 560);
    std::string title = opt.title.empty() ? fp.ligandId : opt.title;
    double y = margin;
    std::ostringstream body;
    if (!title.empty() || !fp.pocketName.empty()) {
        std::string t = title + (fp.pocketName.empty() ? "" : "  |  " + fp.pocketName);
        body << "<text x=\"" << f1(margin) << "\" y=\"" << f1(y + 12) << "\" font-size=\"13\" font-weight=\"bold\" fill=\"#212529\">" << esc(t) << "</text>\n";
        y += 24;
    }
    if (opt.showLegend) {
        double x = margin, rowY = y;
        for (int c = 0; c < FpClassCount; c++) {
            const std::string name = fpClassTitle(c);
            const double w = 16 + name.size() * 5.6 + 14;
            if (x + w > width - margin && x > margin) { x = margin; rowY += 16; }
            body << glyph(c, x + 5, rowY + 6, 3.6) << "<text x=\"" << f1(x + 13) << "\" y=\"" << f1(rowY + 9.5)
                 << "\" font-size=\"9\" fill=\"#495057\">" << esc(name) << "</text>\n";
            x += w;
        }
        y = rowY + 22;
    }

    // ---- rows of regions
    for (const auto& row : rows) {
        double x = margin;
        for (const auto& chunk : row) {
            const double w = chunk.cells.size() * cw;
            body << "<g><title>" << esc(chunk.region) << "</title><rect x=\"" << f1(x) << "\" y=\"" << f1(y) << "\" width=\"" << f1(w - 2)
                 << "\" height=\"" << f1(headerH) << "\" rx=\"2\" fill=\"" << colorOfRegion(chunk.region) << "\"/>\n";
            // the region name must fit in the band: truncate with an ellipsis (the full name is in the tooltip)
            std::string shown = chunk.region;
            const size_t maxChars = static_cast<size_t>(std::max(1.0, (w - 8) / 5.6));
            if (shown.size() > maxChars) shown = maxChars > 1 ? shown.substr(0, maxChars - 1) + "\xE2\x80\xA6" : shown.substr(0, 1);
            body << "<text x=\"" << f1(x + (w - 2) / 2) << "\" y=\"" << f1(y + 11) << "\" font-size=\"9.5\" text-anchor=\"middle\" fill=\"#212529\">"
                 << esc(shown) << "</text></g>\n";
            const double cy0 = y + headerH + 2;
            for (size_t k = 0; k < chunk.cells.size(); k++) {
                const FingerprintEntry& e = *chunk.cells[k];
                const double cx0 = x + k * cw;
                const bool any = e.any();
                std::ostringstream tip;
                tip << e.residue;
                if (!e.position.empty()) tip << " (" << e.position << ")";
                if (any) {
                    tip << ":";
                    for (int c = 0; c < FpClassCount; c++)
                        if (e.counts[c] > 0) tip << " " << fpClassName(c) << " x" << e.counts[c] << " (min " << f1(e.minDistance[c]) << " A);";
                } else tip << ": no interactions";
                body << "<g><title>" << esc(tip.str()) << "</title>\n";
                body << "<rect x=\"" << f1(cx0 + 0.5) << "\" y=\"" << f1(cy0) << "\" width=\"" << f1(cw - 3) << "\" height=\"" << f1(cellH)
                     << "\" rx=\"2\" fill=\"" << (any ? "#f1f3f5" : "#ffffff") << "\" stroke=\"#dee2e6\" stroke-width=\"0.8\"/>\n";
                if (opt.showPositionLabels && !e.position.empty())
                    body << "<text x=\"" << f1(cx0 + (cw - 2) / 2) << "\" y=\"" << f1(cy0 + 8.5) << "\" font-size=\"6.5\" text-anchor=\"middle\" fill=\"#868e96\">"
                         << esc(e.position) << "</text>\n";
                const char letter = e.resName.empty() ? '?' : oneLetterCode(e.resName);
                body << "<text x=\"" << f1(cx0 + (cw - 2) / 2) << "\" y=\"" << f1(cy0 + labelH + 11) << "\" font-size=\"12\" font-weight=\"bold\" text-anchor=\"middle\" fill=\""
                     << letterColor(e) << "\">" << esc(std::string(1, letter)) << "</text>\n";
                body << "<text x=\"" << f1(cx0 + (cw - 2) / 2) << "\" y=\"" << f1(cy0 + labelH + letterH + 8) << "\" font-size=\"7\" text-anchor=\"middle\" fill=\"#868e96\">"
                     << e.resi << (e.icode != ' ' ? std::string(1, e.icode) : "") << "</text>\n";
                int slot = 0;
                const double gy0 = cy0 + labelH + letterH + numberH + 5;
                for (int c = 0; c < FpClassCount; c++) {
                    if (e.counts[c] == 0) continue;
                    body << glyph(c, cx0 + 5.5 + (slot % 4) * 6.3, gy0 + (slot / 4) * 8, 2.6);
                    slot++;
                }
                body << "</g>\n";
            }
            x += w;
        }
        y += headerH + 2 + cellH + gap;
    }
    if (rows.empty()) {
        body << "<text x=\"" << f1(margin) << "\" y=\"" << f1(y + 14) << "\" font-size=\"11\" fill=\"#868e96\">No interactions</text>\n";
        y += 30;
    }

    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << f1(width) << "\" height=\"" << f1(y + margin - gap + 4)
        << "\" viewBox=\"0 0 " << f1(width) << " " << f1(y + margin - gap + 4) << "\" font-family=\"Helvetica, Arial, sans-serif\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
        << body.str() << "</svg>\n";
    return svg.str();
}

}  // namespace ligint
