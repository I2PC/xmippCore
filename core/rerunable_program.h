/***************************************************************************
 *
 * Authors:     David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#pragma once

#include "metadata_db.h"
#include "metadata_vec.h"
#include "xmipp_filename.h"

/**
 * This class adds rerun functionality to Xmipp program.
 * Program is expected to store the progres in the MetaData file stored in the
 *file system. In case of failure or interruption, this file can be loaded and
 *the program will automatically pick up where it left.
 **/
class Rerunable {
protected:
  /**
   * @fn name of the file where intermediary results will be stored
   **/
  Rerunable(const FileName &fn) : fnDone(fn) {}

  /**
   * This method will fill remove processed records from the md.
   * Md will not be changed if no data has been processed.
   **/
  virtual void createWorkFiles(bool resume, MetaData *md) {
    if (nullptr == md) {
      REPORT_ERROR(ERR_MD,
                   "Null pointer passed. "
                   "If you can reproduce this, please contact developers.");
    }
    if (resume && fnDone.exists()) {
      MetaDataDb done(fnDone);
      auto *toDo = dynamic_cast<MetaDataDb *>(md);
      if (nullptr == toDo) {
        MetaDataDb tmp(*md);
        tmp.subtraction(done, MDL_IMAGE);
        *md = tmp;
      } else {
        toDo->subtraction(done, MDL_IMAGE);
      }
    } else // if not exists create metadata only with headers
    {
      MetaDataVec mdDone;
      for (const auto &l : this->getLabelsForEmpty()) {
        mdDone.addLabel(l);
      }
      mdDone.write(fnDone);
    }
  }

  /**
   * Returns labels to be used in the empty progres MetaData file
   **/
  virtual std::vector<MDLabel> getLabelsForEmpty() = 0;

  const FileName &getFileName() const { return this->fnDone; }

  void setFileName(const FileName &fn) { fnDone = fn; }

private:
  FileName fnDone;
};
