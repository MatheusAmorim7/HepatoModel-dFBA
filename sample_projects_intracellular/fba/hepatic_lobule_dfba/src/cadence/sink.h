#include <zlib.h>
#include <streambuf>
struct SinkBuf:std::streambuf {
    std::filebuf file;gzFile gz=nullptr;bool zipped,closed=false;
    SinkBuf(const std::string& path,bool zip):zipped(zip){
        if(zip){gz=gzopen(path.c_str(),"wb1");need(gz!=nullptr,"Sink_open_failed");need(gzbuffer(gz,65536)==0,"Sink_buffer_failed");}
        else need(file.open(path,std::ios::out|std::ios::binary)!=nullptr,"Sink_open_failed");
    }
    std::streamsize xsputn(const char* data,std::streamsize n)override{
        if(closed)return 0;if(!zipped)return file.sputn(data,n);
        std::streamsize done=0;while(done<n){unsigned k=unsigned(std::min<std::streamsize>(n-done,1048576));int wrote=gzwrite(gz,data+done,k);if(wrote!=int(k))return done;done+=wrote;}return done;
    }
    int sync()override{return closed?-1:(zipped?(gzflush(gz,Z_SYNC_FLUSH)==Z_OK?0:-1):file.pubsync());}
    int_type overflow(int_type c)override{if(traits_type::eq_int_type(c,traits_type::eof()))return traits_type::not_eof(c);char x=traits_type::to_char_type(c);return xsputn(&x,1)==1?c:traits_type::eof();}
    void finish(){need(!closed,"Sink_already_closed");closed=true;bool ok=zipped?gzclose(gz)==Z_OK:file.close()!=nullptr;gz=nullptr;need(ok,"Sink_close_failed");}
    ~SinkBuf(){if(!closed){if(zipped&&gz)gzclose(gz);else if(!zipped)file.close();}}
};
struct BinarySink:std::ostream {
    SinkBuf buffer;
    BinarySink(const std::string& p,bool zip):std::ostream(nullptr),buffer(p,zip){rdbuf(&buffer);}
    void finish(){need(bool(*this),"Sink_write_failed");buffer.finish();}
};
int sink_fault_test(const std::string& mode,const std::string& output){
    bool rejected=false;std::string reason;
    try{
        BinarySink s("/dev/full",mode!="raw_write");
        if(mode=="gzip_close"){s.put('x');need(bool(s),"Unexpected_early_close_test_failure");}
        else{std::vector<char> bytes(1048576);unsigned state=1234567;for(char& c:bytes){state^=state<<13;state^=state>>17;state^=state<<5;c=char(state);};s.write(bytes.data(),bytes.size());need(bool(s),"Sink_write_failed");}
        s.finish();
    }catch(const std::exception& e){reason=e.what();rejected=reason==(mode=="gzip_close"?"Sink_close_failed":"Sink_write_failed");}
    need(rejected,"Expected_sink_failure_not_observed");std::ofstream f(output);f<<"{\"expected_failure_detected\":true,\"mode\":\""<<mode<<"\",\"reason\":\""<<reason<<"\",\"simulation_executed\":false}\n";return 0;
}
