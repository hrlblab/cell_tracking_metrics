/*
 * TRAMeasure.cpp
 *
 * Calculate tracking accuracy
 *
 * Martin Maska (xmaska@fi.muni.cz) 2019
 * 
 */

#include <i3d/image3d.h>
#include <i3d/histogram.h>
#include <wx/dir.h>
#include <wx/filename.h>
#include <wx/textfile.h>

/////////////////////////////////////////////////////////////////
//                                                             //
//                         Local Classes                       //
//                                                             //
/////////////////////////////////////////////////////////////////

/**********************************************************************/

/** Track representation. */
template <class T> class Track
{
	public:
		/** Track identificator. */
		T m_id;
		/** The number of frame in which the track begins. */
		size_t m_begin;
		/** The number of frame in which the track ends. */
		size_t m_end;
		/** Identificator of the parent track. */
		T m_parent;
	
	public:
		/** Constructor. */
		Track(T id = 0, size_t begin = -1, size_t end = -1, T parent = 0) : m_id(id), m_begin(begin), m_end(end), m_parent(parent) {};
};

/**********************************************************************/

/** Temporal level representation. */
template <class T> class TemporalLevel
{
	public:
		/** Temporal level. */
		size_t m_level;
		/** List of labels in the reference image. */
		std::vector<std::pair<T, size_t> > m_gt_lab;
		/** List of labels in the computed image. */
		std::vector<std::pair<T, size_t> > m_res_lab;
		/** Matching matrix. */
		std::vector<size_t> m_match;
		/** Indices of reference vertex matching. The value -1 corresponds to an FN vertex. */
		std::vector<size_t> m_gt_match;
		/** Lists of indices of computed vertex matching. The empty list corresponds to an FP vertex. */
		std::vector<std::vector<size_t> > m_res_match;

	public:
		/** Constructor. */
		TemporalLevel(size_t level) : m_level(level), m_gt_lab(), m_res_lab(), m_match(), m_gt_match(), m_res_match() {};
};

/**********************************************************************/

/** Penalty configuration representation. */
class PenaltyConfig
{
	public:
		/** The penalty for a splitting operation. */
		double m_ns;
		/** The penalty for a false negative node. */
		double m_fn;
		/** The penalty for a false positive node. */
		double m_fp;
		/** The penalty for a redundant edge. */
		double m_ed;
		/** The penalty for a missing edge. */
		double m_ea;
		/** The penalty for an edge with wrong semantics. */
		double m_ec;

	public:
		/** Constructor. */
		PenaltyConfig(double ns, double fn, double fp, double ed, double ea, double ec) : m_ns(ns), m_fn(fn), m_fp(fp), m_ed(ed), m_ea(ea), m_ec(ec) 
		{
			if (ns < 0.0 || fn < 0.0 || fp < 0.0 || ed < 0.0 || ea < 0.0 || ec < 0.0)
			{
				printf("All weights must be nonnegative numbers!\n");
				exit(-1);
			}
		};
};

/**********************************************************************/

/////////////////////////////////////////////////////////////////
//                                                             //
//                         Local Functions                     //
//                                                             //
/////////////////////////////////////////////////////////////////

/**********************************************************************/

/** Load image data from a given file. Return false if any error occurs. */
template <class T> static bool LoadImage(i3d::Image3d<T> &img, const char *fname)
{
	try 
	{
		img.ReadImage(fname);
	}
	catch (i3d::IOException &e)
	{
		printf("i3d::IOException: %s\n", e.what.c_str());
		return false;
	}
	catch (i3d::InternalException &e)
	{
		printf("i3d::InternalException: %s\n", e.what.c_str());
		return false;
	}
	catch (std::bad_alloc &)
	{
		printf("Not enough memory\n");
		return false;
	}        
	catch (...)
	{
		printf("Unknown exception caught\n");
		return false;
	}

	return true;
}

/**********************************************************************/

/** Load an image from a given file. Return false if any error occurs. */
template <class T> static bool ReadImage(i3d::Image3d<T> &img, const char *fname)
{
	try 
	{
		i3d::ImgVoxelType vt = i3d::ReadImageType(fname);

		if (vt == i3d::Gray16Voxel)
		{
			return LoadImage(img, fname);
		}
		else
		{
			printf("Unsupported image type (%s)\n", i3d::VoxelTypeToString(vt).c_str());
			return false;
		}
	}
	catch (i3d::IOException &e)
	{
		printf("i3d::IOException: %s\n", e.what.c_str());
		return false;
	}
}

/**********************************************************************/

/** Load a track file from a given file. */
template <class T> static void LoadTrackFile(const wxString &fname, std::map<T, Track<T> > &track_list)
{
	wxTextFile file(fname);
	file.Open();

	if (file.IsOpened())
	{
		wxString line;
		T id, parent;
		size_t begin, end;

		for (size_t i = 0; i < file.GetLineCount(); ++i)
		{
			line = file.GetLine(i);
			std::istringstream str(std::string(line.c_str()));
			str >> id >> begin >> end >> parent;

			if (str.fail() || !str.eof())
			{
				printf("The track file '%s' is in incorrect format! (Line=%d)\n", static_cast<const char *>(fname.c_str()), static_cast<int>(i + 1));
				exit(-1);
			}

			if (track_list.find(id) != track_list.end())
			{
				printf("The track file '%s' is in incorrect format! (Duplicite id=%d)\n", static_cast<const char *>(fname.c_str()), static_cast<int>(id));
				exit(-1);
			}

			if (begin <= end)
			{
				track_list[id] = Track<T>(id, begin, end, parent);
			}
			else
			{
				printf("The track with the label %d is in incorrect format!\n", static_cast<int>(id));
			}
		}
	}
	else
	{
		printf("The track file '%s' cannot be opened!\n", static_cast<const char *>(fname.c_str()));
		exit(-1);
	}
}

/**********************************************************************/

/** Create a list of labels from a given histogram. */
template <class T> static void CreateLabels(const i3d::Histogram &hist, std::vector<std::pair<T, size_t> > &lab, size_t tpos)
{
	lab.clear();

	for (size_t i = 1; i < hist.size(); ++i)
	{
		if (hist[i] > 0)
		{
			lab.push_back(std::pair<T, size_t>(static_cast<T>(i), hist[i]));
		}
	}
}

/**********************************************************************/

/** Find position of a given label. Return -1 if the given label does not exist. */
template <class T> static size_t FindLabel(const std::vector<std::pair<T, size_t> > &lab, T label)
{
	for (size_t i = 0; i < lab.size(); ++i)
	{
		if (lab[i].first == label)
		{
			return i;
		}
	}

	return -1;
}

/**********************************************************************/

/** Go through the given reference and computed images and count their matching grid points. */ 
template <class T> static void CreateMatch(const i3d::Image3d<T> &gt, const i3d::Image3d<T> &res, const std::vector<std::pair<T, size_t> > &gt_lab, const std::vector<std::pair<T, size_t> > &res_lab, std::vector<size_t> &match)
{
	const T *ptr1 = gt.begin();
	const T *ptr2 = res.begin();
	const T *last = gt.end();
	match.resize(gt_lab.size() * res_lab.size(), 0);

	while (ptr1 != last)
	{
		if (*ptr1 > T(0) && *ptr2 > T(0))
		{
			++match[FindLabel(gt_lab, *ptr1) + FindLabel(res_lab, *ptr2) * gt_lab.size()];
		}

		++ptr1;
		++ptr2;
	}
}

/**********************************************************************/

/** Find matching between the reference and computed vertices at a given temporal level. */
template <class T> static void FindMatch(TemporalLevel<T> &level, const PenaltyConfig &penalty, double &aogm, size_t &max_split, std::map<wxString, wxArrayString> &log)
{
	level.m_gt_match.resize(level.m_gt_lab.size(), -1);
	level.m_res_match.resize(level.m_res_lab.size());

	for (size_t i = 0; i < level.m_gt_lab.size(); ++i)
	{
		for (size_t j = 0; j < level.m_res_lab.size(); ++j)
		{
			if (double(level.m_match[j * level.m_gt_lab.size() + i]) / double(level.m_gt_lab[i].second) > 0.5)
			{
				level.m_gt_match[i] = j;
				level.m_res_match[j].push_back(i);
				break;
			}
		}

		if (level.m_gt_match[i] == size_t(-1))
		{
			aogm += penalty.m_fn;
			log["FN"].Add(wxString::Format("T=%d GT_label=%d", static_cast<int>(level.m_level), static_cast<int>(level.m_gt_lab[i].first)));
		}
	}

	size_t num;

	for (size_t i = 0; i < level.m_res_match.size(); ++i)
	{
		num = level.m_res_match[i].size();

		if (num == 0)
		{
			aogm += penalty.m_fp;
			log["FP"].Add(wxString::Format("T=%d Label=%d", static_cast<int>(level.m_level), static_cast<int>(level.m_res_lab[i].first)));
		}
		else if (num > 1)
		{
			aogm += (num - 1) * penalty.m_ns;
			log["NS"].Add(wxString::Format("T=%d Label=%d", static_cast<int>(level.m_level), static_cast<int>(level.m_res_lab[i].first)), num - 1);
			max_split = std::max(max_split, num);
		}
	}
}

/**********************************************************************/

/** Retrieve the temporal position from a given filename and a given number of digits used for encoding the temporal position. */
static bool GetPosition(const wxString &fname, size_t num_digits, size_t &tpos)
{
	if (fname.Length() == (13 + num_digits))
	{
		long tmp;
		
		if (fname.BeforeLast('.').Mid(9).ToCLong(&tmp))
		{
			tpos = static_cast<size_t>(tmp);
			return true;
		}
	}

	return false;
}

/**********************************************************************/

/** Classify labels between the reference and computed images. */
template <class T> static void ClassifyLabels(const char *gt_name, const char *res_name, std::vector<TemporalLevel<T> > &levels, const PenaltyConfig &penalty, double &aogm, size_t &max_split, std::map<wxString, wxArrayString> &log, size_t min_tpos)
{
	// load the reference image
	i3d::Image3d<T> gt_img;
	if (!ReadImage(gt_img, gt_name))
	{
		printf("The image '%s' cannot be opened!\n", gt_name);
		exit(-1);
	}

	// load the computed image
	i3d::Image3d<T> res_img;
	if (!ReadImage(res_img, res_name))
	{
		printf("The image '%s' cannot be opened!\n", res_name);
		exit(-1);
	}
	
	// check the size of both images
	if (gt_img.GetSize() != res_img.GetSize())
	{
		printf("Incompatible image size! (gt='%s', res='%s')\n", gt_name, res_name);
		exit(-1);
	}

	// match labels
	i3d::Histogram gt_hist, res_hist;
	i3d::IntensityHist(gt_img, gt_hist);
	i3d::IntensityHist(res_img, res_hist);

	TemporalLevel<T> level(levels.size() + min_tpos);
	CreateLabels(gt_hist, level.m_gt_lab, level.m_level);
	CreateLabels(res_hist, level.m_res_lab, level.m_level);
	CreateMatch(gt_img, res_img, level.m_gt_lab, level.m_res_lab, level.m_match);
	FindMatch(level, penalty, aogm, max_split, log);
	levels.push_back(level);
}

/**********************************************************************/

/** Find minimum and maximum temporal positions in given tracking results. */
template <class T> static void FindTemporalRange(const std::map<T, Track<T> > &tracks, size_t &min_tpos, size_t &max_tpos)
{
	min_tpos = std::numeric_limits<size_t>::max();
	max_tpos = std::numeric_limits<size_t>::min();

	typename std::map<T, Track<T> >::const_iterator it = tracks.begin();

	while (it != tracks.end())
	{
		min_tpos = std::min(min_tpos, it -> second.m_begin);
		max_tpos = std::max(max_tpos, it -> second.m_end);
		++it;
	}
}

/**********************************************************************/

/** Trim given tracking results to belong to a given temporal range. */
template <class T> static void TrimTracks(std::map<T, Track<T> > &tracks, size_t min_tpos, size_t max_tpos)
{
	// trim parent connections that begin before minimum temporal position
	typename std::map<T, Track<T> >::iterator it = tracks.begin();

	while (it != tracks.end())
	{
		if (it -> second.m_parent > T(0) && tracks.at(it -> second.m_parent).m_end < min_tpos)
		{
			it -> second.m_parent = T(0);
		}

		++it;
	}
	
	// trim tracks to belong to a given temporal range
	typename std::map<T, Track<T> > tmp = tracks;
	tracks.clear();
	it = tmp.begin();

	while (it != tmp.end())
	{
		if (it -> second.m_begin <= max_tpos && it -> second.m_end >= min_tpos)
		{
			tracks[it -> first] = Track<T>(it -> first, std::max(min_tpos, it -> second.m_begin), std::min(max_tpos, it -> second.m_end), it -> second.m_parent);
		}

		++it;
	}
}

/**********************************************************************/

/** Check the consistency of the track files with the image data. */
template <class T> static void CheckConsistency(const std::map<T, Track<T> > &gt_tracks, const std::map<T, Track<T> > &res_tracks, const std::vector<TemporalLevel<T> > &levels, size_t min_tpos)
{
	typename std::map<T, Track<T> >::const_iterator it = gt_tracks.begin();

	while (it != gt_tracks.end())
	{
		if ((it -> second.m_begin - min_tpos) >= levels.size() || (it -> second.m_end - min_tpos) >= levels.size())
		{
			printf("The reference track with label %d is not consistent with the image data!\n", static_cast<int>(it -> first));
			exit(-1);
		}

		if (it -> second.m_parent > T(0) && gt_tracks.find(it -> second.m_parent) == gt_tracks.end())
		{
			printf("Non-existing parent track with label %d for the reference track with label %d!\n", static_cast<int>(it -> second.m_parent), static_cast<int>(it -> first));
			exit(-1);
		}

		if (it -> second.m_parent > T(0) && it -> second.m_begin <= gt_tracks.at(it -> second.m_parent).m_end)
		{
			printf("Invalid parent connection for the reference track with label %d!\n", static_cast<int>(it -> first));
			exit(-1);
		}

		for (size_t i = it -> second.m_begin; i <= it -> second.m_end; ++i)
		{
			if (FindLabel(levels[i - min_tpos].m_gt_lab, it -> first) == size_t(-1))
			{
				printf("The reference track with label %d is not consistent with the image data!\n", static_cast<int>(it -> first));
				exit(-1);
			}
		}

		++it;
	}

	it = res_tracks.begin();

	while (it != res_tracks.end())
	{
		if ((it -> second.m_begin - min_tpos) >= levels.size() || (it -> second.m_end - min_tpos) >= levels.size())
		{
			printf("The computed track with label %d is not consistent with the image data!\n", static_cast<int>(it -> first));
			exit(-1);
		}

		if (it -> second.m_parent > T(0) && res_tracks.find(it -> second.m_parent) == res_tracks.end())
		{
			printf("Non-existing parent track with label %d for the computed track with label %d!\n", static_cast<int>(it -> second.m_parent), static_cast<int>(it -> first));
			exit(-1);
		}

		if (it -> second.m_parent > T(0) && it -> second.m_begin <= res_tracks.at(it -> second.m_parent).m_end)
		{
			printf("Invalid parent connection for the computed track with label %d!\n", static_cast<int>(it -> first));
			exit(-1);
		}

		for (size_t i = it -> second.m_begin; i <= it -> second.m_end; ++i)
		{
			if (FindLabel(levels[i - min_tpos].m_res_lab, it -> first) == size_t(-1))
			{
				printf("The computed track with label %d is not consistent with the image data!\n", static_cast<int>(it -> first));
				exit(-1);
			}
		}

		++it;
	}

	for (size_t i = 0; i < levels.size(); ++i)
	{
		for (size_t j = 0; j < levels[i].m_gt_lab.size(); ++j)
		{
			it = gt_tracks.find(levels[i].m_gt_lab[j].first);

			if (it == gt_tracks.end() || i < (it -> second.m_begin - min_tpos) || i > (it -> second.m_end - min_tpos))
			{
				printf("The reference track with label %d is not consistent with the image data!\n", static_cast<int>(levels[i].m_gt_lab[j].first));
				exit(-1);
			}
		}

		for (size_t j = 0; j < levels[i].m_res_lab.size(); ++j)
		{
			it = res_tracks.find(levels[i].m_res_lab[j].first);

			if (it == res_tracks.end() || i < (it -> second.m_begin - min_tpos) || i > (it -> second.m_end - min_tpos))
			{
				printf("The computed track with label %d is not consistent with the image data!\n", static_cast<int>(levels[i].m_res_lab[j].first));
				exit(-1);
			}
		}
	}
}

/**********************************************************************/

/** Get an index of a reference matching vertex for a given label at a given temporal level. */
template <class T> static size_t GetGTMatch(const TemporalLevel<T> &level, T label)
{
	return level.m_gt_match[FindLabel(level.m_gt_lab, label)];
}

/**********************************************************************/

/** Get an index of a computed matching vertex for a given label at a given temporal level. */
template <class T> static std::vector<size_t> GetResMatch(const TemporalLevel<T> &level, T label)
{
	size_t index = FindLabel(level.m_res_lab, label);

	if (index == size_t(-1))
	{
		std::vector<size_t> tmp(1, size_t(-1));
		return tmp;
	}
	else
	{
		return level.m_res_match[index];
	}
}

/**********************************************************************/

/** Check if there is an edge of a given type between given temporal levels in the reference tracks. */
template <class T> static bool ExistGTEdge(const std::vector<TemporalLevel<T> > &levels, size_t start_level, size_t start_index, size_t end_level, size_t end_index, const std::map<T, Track<T> > &tracks, bool &parent, size_t min_tpos)
{
	if (start_index != size_t(-1) && end_index != size_t(-1))
	{
		// check the existence and type of the edge
		T start_label = levels[start_level].m_gt_lab[start_index].first;
		T end_label = levels[end_level].m_gt_lab[end_index].first;

		if (start_label == end_label)
		{
			// the edge is within a track
			if ((start_level + 1) == end_level)
			{
				parent = false;
				return true;
			}
		}
		else
		{
			// the edge connects two tracks
			const Track<T> &parent_node = tracks.at(start_label);
			const Track<T> &child_node = tracks.at(end_label);

			if ((parent_node.m_end - min_tpos) == start_level && (child_node.m_begin - min_tpos) == end_level && child_node.m_parent == start_label)
			{
				parent = true;
				return true;
			}
		}
	}

	return false;
}

/**********************************************************************/

/** Check if there is an edge between given temporal levels in the computed tracks. */
template <class T> static bool ExistResEdge(const std::vector<TemporalLevel<T> > &levels, size_t start_level, size_t start_index, size_t end_level, size_t end_index, const std::map<T, Track<T> > &tracks, size_t min_tpos)
{
	if (start_index != size_t(-1) && end_index != size_t(-1) && levels[start_level].m_res_match[start_index].size() == 1 && levels[end_level].m_res_match[end_index].size() == 1)
	{
		T start_label = levels[start_level].m_res_lab[start_index].first;
		T end_label = levels[end_level].m_res_lab[end_index].first;
	
		if (start_label == end_label)
		{
			// the edge is within a track
			return ((start_level + 1) == end_level);
		}
		else
		{
			// the edge connects two tracks
			const Track<T> &parent = tracks.at(start_label);
			const Track<T> &child = tracks.at(end_label);
	
			return ((parent.m_end - min_tpos) == start_level && (child.m_begin - min_tpos) == end_level && child.m_parent == start_label);
		}
	}

	return false;
}

/**********************************************************************/

/** Find edges in the computed tracks that must be removed or altered. */
template <class T> static void FindEDAndECEdges(const std::vector<TemporalLevel<T> > &levels, const std::map<T, Track<T> > &gt_tracks, const std::map<T, Track<T> > &res_tracks, const PenaltyConfig &penalty, double &aogm, std::map<wxString, wxArrayString> &log, size_t min_tpos)
{
	bool parent;
	size_t start_level, end_level;
	std::vector<size_t> start_match, end_match;
	typename std::map<T, Track<T> >::const_iterator it = res_tracks.begin();

	while (it != res_tracks.end())
	{
		// check the edge between the first node of the current track and the last one of the parent track
		end_level = it -> second.m_begin - min_tpos;
		end_match = GetResMatch(levels[end_level], it -> first);

		if (it -> second.m_parent > T(0))
		{
			start_level = res_tracks.at(it -> second.m_parent).m_end - min_tpos;
			start_match = GetResMatch(levels[start_level], it -> second.m_parent);

			if (start_match.size() == 1 && end_match.size() == 1)
			{
				if (ExistGTEdge(levels, start_level, start_match[0], end_level, end_match[0], gt_tracks, parent, min_tpos))
				{
					if (!parent)
					{
						aogm += penalty.m_ec;
						log["EC"].Add(wxString::Format("[T=%d Label=%d] -> [T=%d Label=%d]", static_cast<int>(start_level + min_tpos), static_cast<int>(it->second.m_parent), static_cast<int>(end_level + min_tpos), static_cast<int>(it->first)));
					}
				}
				else
				{
					aogm += penalty.m_ed;
					log["ED"].Add(wxString::Format("[T=%d Label=%d] -> [T=%d Label=%d]", static_cast<int>(start_level + min_tpos), static_cast<int>(it->second.m_parent), static_cast<int>(end_level + min_tpos), static_cast<int>(it->first)));
				}
			}
		}

		// check edges within the current track
		for (size_t i = it -> second.m_begin; i < it -> second.m_end; ++i)
		{
			start_level = end_level;
			start_match = end_match;
			end_level = i + 1 - min_tpos;
			end_match = GetResMatch(levels[end_level], it -> first);

			if (start_match.size() == 1 && end_match.size() == 1)
			{
				if (ExistGTEdge(levels, start_level, start_match[0], end_level, end_match[0], gt_tracks, parent, min_tpos))
				{
					if (parent)
					{
						aogm += penalty.m_ec;
						log["EC"].Add(wxString::Format("[T=%d Label=%d] -> [T=%d Label=%d]", static_cast<int>(start_level + min_tpos), static_cast<int>(it->first), static_cast<int>(end_level + min_tpos), static_cast<int>(it->first)));
					}
				}
				else
				{
					aogm += penalty.m_ed;
					log["ED"].Add(wxString::Format("[T=%d Label=%d] -> [T=%d Label=%d]", static_cast<int>(start_level + min_tpos), static_cast<int>(it->first), static_cast<int>(end_level + min_tpos), static_cast<int>(it->first)));
				}
			}
		}

		++it;
	}
}

/**********************************************************************/

/** Find edges in the reference tracks that must be added. */
template <class T> static void FindEAEdges(const std::vector<TemporalLevel<T> > &levels, const std::map<T, Track<T> > &gt_tracks, const std::map<T, Track<T> > &res_tracks, const PenaltyConfig &penalty, double &aogm, std::map<wxString, wxArrayString> &log, size_t min_tpos)
{
	size_t start_level, end_level, start_index, end_index;
	typename std::map<T, Track<T> >::const_iterator it = gt_tracks.begin();

	while (it != gt_tracks.end())
	{
		// check the edge between the first node of the current track and the last one of the parent track
		end_level = it -> second.m_begin - min_tpos;
		end_index = GetGTMatch(levels[end_level], it -> first);

		if (it -> second.m_parent > T(0))
		{
			start_level = gt_tracks.at(it -> second.m_parent).m_end - min_tpos;
			start_index = GetGTMatch(levels[start_level], it -> second.m_parent);

			if (!ExistResEdge(levels, start_level, start_index, end_level, end_index, res_tracks, min_tpos))
			{
				aogm += penalty.m_ea;
				log["EA"].Add(wxString::Format("[T=%d GT_label=%d] -> [T=%d GT_label=%d]", static_cast<int>(start_level + min_tpos), static_cast<int>(it -> second.m_parent), static_cast<int>(end_level + min_tpos), static_cast<int>(it -> first)));
			}
		}

		// check edges within the current track
		for (size_t i = it -> second.m_begin; i < it -> second.m_end; ++i)
		{
			start_level = end_level;
			start_index = end_index;
			end_level = i + 1 - min_tpos;
			end_index = GetGTMatch(levels[end_level], it -> first);

			if (!ExistResEdge(levels, start_level, start_index, end_level, end_index, res_tracks, min_tpos))
			{
				aogm += penalty.m_ea;
				log["EA"].Add(wxString::Format("[T=%d GT_label=%d] -> [T=%d GT_label=%d]", static_cast<int>(start_level + min_tpos), static_cast<int>(it -> first), static_cast<int>(end_level + min_tpos), static_cast<int>(it -> first)));
			}
		}

		++it;
	}
}

/**********************************************************************/

/** Log errors of a given class to a given file. */
static void LogClass(const wxString &class_name, wxFile &file, std::map<wxString, wxArrayString> &log)
{
	const wxArrayString &list = log[class_name];
	
	for(size_t i = 0; i < list.size(); ++i)
	{
		file.Write(list[i] + wxTextFile::GetEOL());
	}
}

/**********************************************************************/

/////////////////////////////////////////////////////////////////
//                                                             //
//                           Entry Point                       //
//                                                             //
/////////////////////////////////////////////////////////////////

/**********************************************************************/

/** Entry point of the console application. */
int main(int argc, char **argv)
{
	if (argc != 4)
	{
		 std::cout << "A simple console application intended for calculating the tracking accuracy\n\n"
	                  "Copyright (C) 2019 Centre for Biomedical Image Analysis, Masaryk University\n"
					  "This is free software; see the source for copying conditions. There is NO\nwarranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n";
    
		 std::cout << "Usage: "<< argv[0] << " dir seq num_digits\n";

		 return -1;
	}
	
	wxString gt_pref = wxString::Format("%s/%s_GT/TRA/", argv[1], argv[2]);
	wxString res_pref = wxString::Format("%s/%s_RES/", argv[1], argv[2]);
		
	if (wxDir::Exists(gt_pref))
	{
		if (wxDir::Exists(res_pref))
		{
			wxString res_fname_temp;
			int num_digits = atoi(argv[3]);

			if (num_digits > 2)
			{
				res_fname_temp = wxString::Format("%smask%%0%dd.tif", res_pref, num_digits);
			}
			else
			{
				printf("The number of digits used for encoding temporal positions must be at least 3! (num_digits=%d)\n", num_digits);
				return -1;
			}

			PenaltyConfig penalty(5.0, 10.0, 1.0, 1.0, 1.5, 1.0);
			double aogm = 0.0;

			std::map<i3d::GRAY16, Track<i3d::GRAY16> > gt_tracks, res_tracks;
			LoadTrackFile(gt_pref + "man_track.txt", gt_tracks);
			LoadTrackFile(res_pref + "res_track.txt", res_tracks);
			
			size_t min_tpos, max_tpos;
			FindTemporalRange(gt_tracks, min_tpos, max_tpos);
			
			std::vector<TemporalLevel<i3d::GRAY16> > levels;
			wxArrayString files;
			std::map<wxString, wxArrayString> log;
			log["NS"].Add(wxString::Format("----------Splitting Operations (Penalty=%g)----------", penalty.m_ns));
			log["FN"].Add(wxString::Format("----------False Negative Vertices (Penalty=%g)----------", penalty.m_fn));
			log["FP"].Add(wxString::Format("----------False Positive Vertices (Penalty=%g)----------", penalty.m_fp));
			log["ED"].Add(wxString::Format("----------Redundant Edges To Be Deleted (Penalty=%g)----------", penalty.m_ed));
			log["EA"].Add(wxString::Format("----------Edges To Be Added (Penalty=%g)----------", penalty.m_ea));
			log["EC"].Add(wxString::Format("----------Edges with Wrong Semantics (Penalty=%g)----------", penalty.m_ec));
			size_t tpos = -1;
			size_t max_split = 1;

			wxDir::GetAllFiles(gt_pref, &files, "man_track*.tif", wxDIR_FILES);

			// sort the listed files alphabetically
			files.Sort();

			for (size_t i = 0; i < files.size(); ++i)
			{
				if (GetPosition(wxFileName(files[i]).GetFullName(), num_digits, tpos))
				{
					ClassifyLabels(static_cast<const char *>(files[i].c_str()), static_cast<const char *>(wxString::Format(res_fname_temp, static_cast<int>(tpos)).c_str()), levels, penalty, aogm, max_split, log, min_tpos);
				}
			}

			if (levels.size() == 0)
			{
				printf("No reference image found!\n");
				return -1;
			}
			else
			{
				CheckConsistency(gt_tracks, res_tracks, levels, min_tpos);
				TrimTracks(res_tracks, min_tpos, max_tpos);

				if (gt_tracks.size() == 0)
				{
					printf("No reference marker found!\n");
					return -1;
				}

				// check the minimality condition
				if ((max_split - 1) * penalty.m_ns > (penalty.m_fp + max_split * penalty.m_fn))
				{
					printf("Warning: The minimality condition broken! (m*=%d)\n", static_cast<int>(max_split));
				}
				
				FindEDAndECEdges(levels, gt_tracks, res_tracks, penalty, aogm, log, min_tpos);
				FindEAEdges(levels, gt_tracks, res_tracks, penalty, aogm, log, min_tpos);

				size_t num_par = 0, sum = 0;
				std::map<i3d::GRAY16, Track<i3d::GRAY16> >::const_iterator it = gt_tracks.begin();

				while (it != gt_tracks.end())
				{
					sum += it -> second.m_end - it -> second.m_begin;

					if (it -> second.m_parent > i3d::GRAY16(0))
					{
						++num_par;
					}

					++it;
				}

				// display the final value
				double val_empty = penalty.m_fn * (sum + gt_tracks.size()) + penalty.m_ea * (sum + num_par); 
				double val = 1.0 - std::min(aogm, val_empty) / val_empty;

				// display the TRA value
				printf("TRA measure: %.6f\n", val);

				// create the log file
				wxFile log_file;
				wxString log_name = res_pref + "TRA_log.txt";

				if (log_file.Create(log_name, true))
				{
					LogClass("NS", log_file, log);
					LogClass("FN", log_file, log);
					LogClass("FP", log_file, log);
					LogClass("ED", log_file, log);
					LogClass("EA", log_file, log);
					LogClass("EC", log_file, log);

					// log the AOGM value
					log_file.Write(wxString("=================================================================================") + wxTextFile::GetEOL());

					if ((max_split - 1) * penalty.m_ns > (penalty.m_fp + max_split * penalty.m_fn))
					{
						log_file.Write(wxString::Format("Warning: The minimality condition broken! (m*=%d)", static_cast<int>(max_split)) + wxTextFile::GetEOL());
					}
					
					log_file.Write(wxString::Format("TRA measure: %.6f", val) + wxTextFile::GetEOL());
				}
				else
				{
					printf("The log file '%s' cannot be created!\n", static_cast<const char *>(log_name.c_str()));
					return -1;
				}
			}
		}
		else
		{
			printf("The directory '%s' does not exist!\n", static_cast<const char *>(res_pref.c_str()));
			return -1;
		}
	}
	else
	{
		printf("The directory '%s' does not exist!\n", static_cast<const char *>(gt_pref.c_str()));
		return -1;
	}
			
	return 0;
}

