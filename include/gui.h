#pragma once

#include <QFile>
#include <QFileDialog>
#include <QImage>
#include <QLabel>
#include <QMainWindow>
#include <QPushButton>
#include <QToolBar>

#include "image.h"
#include "raytracer.h"
#include "viewer.h"

class Gui : public QMainWindow {
  public:
    Gui() = delete;

    Gui(int width, int height, RayTracer raytracer, QWindow* parent = nullptr) {

        QToolBar* toolbar = new QToolBar(this);
        toolbar->setMovable(false);

        _saveButton = new QPushButton("Save as ...", this);
        toolbar->addWidget(_saveButton);

        QWidget* spacer = new QWidget();
        spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        toolbar->addWidget(spacer);

        _durationText = new QLabel(this);
        toolbar->addWidget(_durationText);
        this->addToolBar(toolbar);

        _viewer = new Viewer(raytracer, _durationText, this);
        _viewer->resize(width, height);
        this->setCentralWidget(_viewer);

        connect(_saveButton, &QPushButton::clicked, [this]() {
            QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), "render.png",
                                                            tr("Images (*.png);;All Files (*)"));
            QFile file(filename);
            file.open(QIODevice::WriteOnly);
            _viewer->getImage().save(&file, "PNG");
        });

        this->resize(width, height);
    }

    ~Gui() {
        _viewer->stop_raytrace();
    }

  private:
    QLabel* _durationText;
    QPushButton* _saveButton;
    Viewer* _viewer;
};
