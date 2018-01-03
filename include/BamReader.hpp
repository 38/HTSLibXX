#ifndef __HTSLIBPP_BAMREADER_HPP__
#define __HTSLIBPP_BAMREADER_HPP__
#include <BamAlignment.hpp>
#include <SamHeader.hpp>
#include <htslib/sam.h>
#include <stdint.h>
#include <string>
#include <queue>

namespace BamTools {
	class BamReader {
		struct _SamFile {
			_SamFile(samFile* fp, uint32_t _idx, hts_idx_t* _ip = NULL) : fp(fp), idx(_idx), ip(_ip) {}
			void destory() 
			{
				if(nullptr != ip) hts_idx_destroy(ip);
				if(nullptr != fp) sam_close(fp);
			}
			bool load_index(const char* filename) 
			{
				if(ip == NULL && NULL == (ip = sam_index_load(fp, filename)))
					return false;
				return true;
			}
			samFile* fp;
			uint32_t idx;
			hts_idx_t* ip;
		};

		struct _MetaData {
			_SamFile file;
			uint32_t size;
			_MetaData(_SamFile _file, uint32_t size) : file(_file.fp, _file.idx, _file.ip)
			{
				this->size = size;
			}
		};

		struct _Comp {
			bool operator()(const std::pair<_MetaData, bam1_t*>& left, const std::pair<_MetaData, bam1_t*>& right) 
			{
				if(left.second->core.tid == right.second->core.tid)
				{
					return left.second->core.pos >= right.second->core.pos;
				}
				return left.second->core.tid >= right.second->core.tid;
			}
		};

		std::vector<_SamFile> _files;
		std::vector<SamHeader> _hdrs;
		std::priority_queue<std::pair<_MetaData, bam1_t*>, std::vector<std::pair<_MetaData, bam1_t*> >, _Comp> _queue;

		bool _read_sam_file(_SamFile file)
		{
			bam1_t *rec_ptr = bam_init1();
			int read_rc;
			if ((read_rc = sam_read1(file.fp, _hdrs[file.idx].GetHeaderStruct(), rec_ptr)) >= 0)
			{
				_queue.push(std::make_pair(_MetaData(file, (uint32_t)(read_rc - 4)) , rec_ptr));
				return true;
			}
			bam_destroy1(rec_ptr);
			return false;
		}

	public:
		bool Open(const std::vector<std::string>& filenames) 
		{
			uint32_t i = 0;
			for(const auto& filename : filenames)
			{
				
				samFile* fp = sam_open(filename.empty() || filename == "stdin" ? "-" : filename.c_str(), "rb");
				if(nullptr == fp) return false;
				bam_hdr_t* hdr = sam_hdr_read(fp);
				if(nullptr == hdr)
					return false;

				_files.push_back(_SamFile(fp, i ++, NULL));
				_hdrs.push_back(SamHeader(filename, hdr));

				_read_sam_file(_files[_files.size() - 1]);
			}

			return true;
		}

		bool Open(const std::string& filename) 
		{
			std::vector<std::string> vec;
			vec.push_back(filename);
			return Open(vec);
		}

		bool IsOpen() const 
		{
			return true;
		}

		const RefVector GetReferenceData() const
		{
			return _hdrs.at(0).GetReferenceData();
		}

		SamHeader GetHeader(int idx = 0) const
		{
			return _hdrs.at(idx);
		}

		std::string GetHeaderText(int idx = -1) const 
		{
			return _hdrs.at(idx).GetHeaderText();
		}
		std::string GetErrorString() const
		{
			return "FIXME: error string is not supported";
		}

		bool GetNextAlignment(BamAlignment& alignment)
		{
			if(_queue.empty()) return false;

			auto& top = _queue.top();

			_SamFile fp = top.first.file;
			
			alignment(_hdrs[fp.idx].Filename(), top.second, top.first.size);

			_queue.pop();

			_read_sam_file(fp);

			return true;
		}

		bool GetNextAlignmentCore(BamAlignment& alignment)
		{
			return GetNextAlignment(alignment);
		}

		void Close(void)
		{
			for(auto& sam : _files)
			{
				sam.destory();
			}

			for(auto& hdr : _hdrs)
			{
				hdr.destory();
			}

			while(!_queue.empty())
			{
				auto& item = _queue.top();
				bam_destroy1(item.second);
				_queue.pop();
			}
		}

		int GetReferenceID(const std::string& refname)
		{
			if(_hdrs.size() == 0) return -1;

			bam_hdr_t* bh = _hdrs[0].GetHeaderStruct();

			for(int i = 0; i < bh->n_targets; i ++)
				if(strncmp(refname.c_str(), bh->target_name[i], bh->target_len[i]) && refname.c_str()[bh->target_len[i]] == 0)
					return i;
			return -1;
		}

		void LocateIndexes() 
		{
			for(auto& sam: _files)
			{
				if(!sam.load_index(_hdrs[sam.idx].Filename()))
				{
					/* TODO(haohou): Load failure */
				}
			}
		}

		bool HasIndexes()
		{
			for(auto& sam: _files)
			{
				if(NULL == sam.ip)
					return false;
			}
			return true;
		}

		bool SetRegion(BamRegion& region)
		{
			/* TODO(haohou): Implemnet all the index related functions */
			return false;
		}


	};
}
#endif
