#pragma once

#include <QImage>

#include <glm/glm.hpp>

struct Image {
    Image() = delete;
    Image(int width, int height) : _image(width, height, QImage::Format_RGB888) { this->clear(); }

    int width() const { return _image.width(); }
    int height() const { return _image.height(); }

    void setPixel(int x, int y, glm::dvec3 c) {
        _image.setPixel(x, y, QColor((int)(255 * c.r), (int)(255 * c.g), (int)(255 * c.b)).rgb());
    }

    glm::dvec3 getPixel(int x, int y) const {
        auto p = _image.pixel(x, y);
        return {qRed(p) / 255., qGreen(p) / 255., qBlue(p) / 255.};
    }

    void clear() { _image.fill(Qt::black); };

  private:
    QImage _image;

    friend class Viewer;
};
