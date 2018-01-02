
struct _Bin_t {
    _Bin_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _Bin_t& that) {_ptr = that._ptr;}
    operator uint32_t() const {return (uint32_t)(_ptr == NULL?0:_ptr->core.bin);}
    const uint32_t& operator=(const uint32_t& val) {
        if(NULL != _ptr) _ptr->core.bin = (uint32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} Bin;

struct _MapQuality_t {
    _MapQuality_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _MapQuality_t& that) {_ptr = that._ptr;}
    operator int16_t() const {return (int16_t)(_ptr == NULL?0:_ptr->core.qual);}
    const int16_t& operator=(const int16_t& val) {
        if(NULL != _ptr) _ptr->core.qual = (int16_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} MapQuality;

struct _Length_t {
    _Length_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _Length_t& that) {_ptr = that._ptr;}
    operator int32_t() const {return (int32_t)(_ptr == NULL?0:_ptr->core.l_qseq);}
    const int32_t& operator=(const int32_t& val) {
        if(NULL != _ptr) _ptr->core.l_qseq = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} Length;

struct _InsertSize_t {
    _InsertSize_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _InsertSize_t& that) {_ptr = that._ptr;}
    operator int32_t() const {return (int32_t)(_ptr == NULL?0:_ptr->core.isize);}
    const int32_t& operator=(const int32_t& val) {
        if(NULL != _ptr) _ptr->core.isize = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} InsertSize;

struct _MatePosition_t {
    _MatePosition_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _MatePosition_t& that) {_ptr = that._ptr;}
    operator int32_t() const {return (int32_t)(_ptr == NULL?0:_ptr->core.mpos);}
    const int32_t& operator=(const int32_t& val) {
        if(NULL != _ptr) _ptr->core.mpos = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} MatePosition;

struct _Position_t {
    _Position_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _Position_t& that) {_ptr = that._ptr;}
    operator int32_t() const {return (int32_t)(_ptr == NULL?0:_ptr->core.pos);}
    const int32_t& operator=(const int32_t& val) {
        if(NULL != _ptr) _ptr->core.pos = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} Position;

struct _MateRefID_t {
    _MateRefID_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _MateRefID_t& that) {_ptr = that._ptr;}
    operator int32_t() const {return (int32_t)(_ptr == NULL?0:_ptr->core.mtid);}
    const int32_t& operator=(const int32_t& val) {
        if(NULL != _ptr) _ptr->core.mtid = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} MateRefID;

struct _RefID_t {
    _RefID_t() : _ptr(NULL) {}
    void set(bam1_t* ptr) {_ptr = ptr;}
    void set(const _RefID_t& that) {_ptr = that._ptr;}
    operator int32_t() const {return (int32_t)(_ptr == NULL?0:_ptr->core.tid);}
    const int32_t& operator=(const int32_t& val) {
        if(NULL != _ptr) _ptr->core.tid = (int32_t)val;
        return val;
    }
    private:
        bam1_t* _ptr;
} RefID;
void setup(bam1_t* bam)
{
    Bin.set(bam);
    MapQuality.set(bam);
    Length.set(bam);
    InsertSize.set(bam);
    MatePosition.set(bam);
    Position.set(bam);
    MateRefID.set(bam);
    RefID.set(bam);
}
void setup(const BamAlignment& bam)
{
    Bin.set(bam.Bin);
    MapQuality.set(bam.MapQuality);
    Length.set(bam.Length);
    InsertSize.set(bam.InsertSize);
    MatePosition.set(bam.MatePosition);
    Position.set(bam.Position);
    MateRefID.set(bam.MateRefID);
    RefID.set(bam.RefID);
}
