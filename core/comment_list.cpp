/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
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

#include "comment_list.h"

void CommentList::addComment(const String &comment, int visible, bool verbatim)
{
    comments.push_back(comment);
    visibility.push_back(visible);
    wikiVerbatim.push_back(verbatim);
}
void CommentList::addComment(const char * comment, bool verbatim)
{
    size_t t=0;
    while(comment[t]=='+' && comment[t]!='\0')
        t++;
    addComment(comment+t,t,verbatim);
}

void CommentList::clear()
{
    comments.clear();
    visibility.clear();
    wikiVerbatim.clear();
}
size_t CommentList::size() const
{
    return comments.size();
}
