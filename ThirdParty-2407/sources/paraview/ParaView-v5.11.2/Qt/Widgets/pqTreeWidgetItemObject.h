/*=========================================================================

   Program: ParaView
   Module:    pqTreeWidgetItemObject.h

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef pqTreeWidgetItemObject_h
#define pqTreeWidgetItemObject_h

#include "pqWidgetsModule.h"
#include <QObject>
#include <QTreeWidgetItem>

/**
 * QTreeWidgetItem subclass with additional signals, slots, and properties
 */
class PQWIDGETS_EXPORT pqTreeWidgetItemObject
  : public QObject
  , public QTreeWidgetItem
{
  Q_OBJECT
  Q_PROPERTY(bool checked READ isChecked WRITE setChecked)
public:
  /**
   * construct list widget item to for QTreeWidget with a string
   */
  pqTreeWidgetItemObject(const QStringList& t, int type = QTreeWidgetItem::UserType);
  pqTreeWidgetItemObject(
    QTreeWidget* p, const QStringList& t, int type = QTreeWidgetItem::UserType);
  pqTreeWidgetItemObject(
    QTreeWidgetItem* p, const QStringList& t, int type = QTreeWidgetItem::UserType);

  /**
   * overload setData() to emit changed signal
   */
  void setData(int column, int role, const QVariant& v) override;

public Q_SLOTS: // NOLINT(readability-redundant-access-specifiers)
  /**
   * get the check true/false
   */
  bool isChecked() const;
  /**
   * set the check state true/false
   */
  void setChecked(bool v);

Q_SIGNALS:
  /**
   * signal check state changed
   */
  void checkedStateChanged(bool);

  /**
   * Fired every time setData is called.
   */
  void modified();
};

#endif