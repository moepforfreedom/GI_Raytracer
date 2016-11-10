#pragma once

#include <chrono>
#include <thread>

#include <QPainter>
#include <QTimer>
#include <QWidget>
#include <QCloseEvent>

#include "image.h"
#include "raytracer.h"

class Viewer : public QWidget {
  public:
    Viewer(RayTracer raytracer, QLabel* durationText, QWidget* parent)
        : QWidget(parent), _durationText(durationText), _raytracer(raytracer) {
        _timer = new QTimer(this);
        _timer->setInterval(32);
        _timer->start();
        connect(_timer, &QTimer::timeout, [&]() { this->repaint(); });
        restart_raytrace();
    }

    ~Viewer() {
        stop_raytrace();
    }

    void stop_raytrace() {
        if (_raytracer.running()) {
            _raytracer.stop();
            _thread.join();
        }
    }

    void paintEvent(QPaintEvent*) {
        QPainter painter(this);
        painter.drawImage(QPoint(0, 0), getImage());
    }

    void resizeEvent(QResizeEvent*) { restart_raytrace(); }

    QImage getImage() const { return _raytracer.getImage()->_image; }

  private:
    void restart_raytrace() {
        if (_raytracer.running()) {
            _raytracer.stop();
            _thread.join();
        }
        _raytracer.start();
        _thread = std::thread([&]() {
            _durationText->setText("Running...");
            using namespace std::chrono;
            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            this->_raytracer.run(this->width(), this->height());
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(t2 - t1).count();
            _durationText->setText(QString::number(duration / (double)1000) + " seconds");
        });
    }

    QTimer* _timer;
    QLabel* _durationText;
    RayTracer _raytracer;
    std::thread _thread;
};
