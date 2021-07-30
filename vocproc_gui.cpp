/* 
    LV2 plugin for pitch shifting, pitch correction, vocoding and harmonizing
    singing voice

    (c) Igor Brkic 2010

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


#include <gtkmm.h>
#include <lv2gui.hpp>
#include "vocproc.peg"

#include <iostream>

#include <string.h>


using namespace sigc;
using namespace Gtk;
using namespace std;


// maybe a bit messy but it works
class ScalesHelper {
public:
    ScalesHelper(){}
    ~ScalesHelper(){}
protected:
    void scale_change();
};


class VocProcGUI : public LV2::GUI<VocProcGUI>, public ScalesHelper {
public:

    VocProcGUI(const std::string& URI) {

        refBuilder = Builder::create();

        try
        {
            char tmp2[1024];

            strcpy(tmp2, bundle_path());
            strcat(tmp2, "vocproc_gui.ui");
            refBuilder->add_from_file(tmp2);
        }
        catch(const Glib::FileError& ex)
        {
            cerr << "FileError: " << ex.what() << endl;
        }
        catch(const BuilderError& ex)
        {
            cerr << "BuilderError: " << ex.what() << endl;
        }

        HBox *main_box;
        refBuilder->get_widget("main_box", main_box);
        if(main_box)
        {
            
            // -----------------

            refBuilder->get_widget("sPitchFactor", s_pitch_factor);
            s_pitch_factor->set_range(-12, 12);
            s_pitch_factor->set_increments(1, 2);
            slot<void> slot_pitch = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_pitch_factor),
                mem_fun(*s_pitch_factor, &VScale::get_value));
            s_pitch_factor->signal_value_changed().connect(slot_pitch);

            refBuilder->get_widget("sEffect", s_effect);
            s_effect->set_range(0, 1);
            s_effect->set_increments(0.01, 0.2);
            slot<void> slot_effect = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_effect),
                mem_fun(*s_effect, &VScale::get_value));
            s_effect->signal_value_changed().connect(slot_effect);

            // -----------------

            refBuilder->get_widget("cForVoc", c_for_voc);
            slot<void> slot_for_voc_switch = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_fc_voc_switch),
                mem_fun(*c_for_voc, &CheckButton::get_active));
            c_for_voc->signal_toggled().connect(slot_for_voc_switch);

            HBox *b;
#ifndef NO_VOCODER
            refBuilder->get_widget("hbox3", b);
            d_for_voc=new ComboBoxText();
            d_for_voc->append_text("Formant correction");
            d_for_voc->append_text("Vocoder");
            d_for_voc->set_active(0);
            b->pack_start(*d_for_voc);
            d_for_voc->show();
            slot<void> slot_for_voc = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_fc_voc),
                mem_fun(*d_for_voc, &ComboBoxText::get_active_row_number));
            d_for_voc->signal_changed().connect(slot_for_voc);
#else
            Label *tmp_l;
            refBuilder->get_widget("label2", tmp_l);
            tmp_l->set_markup("<b>formant correction</b>");
#endif

            // -----------------
            
            refBuilder->get_widget("cAutoTune", c_autotune);
            slot<void> slot_autotune = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_pitch_correction),
                mem_fun(*c_autotune, &CheckButton::get_active));
            c_autotune->signal_toggled().connect(slot_autotune);

            // scales
            refBuilder->get_widget("hbox4", b);
            d_scales_key=new ComboBoxText();
            d_scales_mode=new ComboBoxText();


            d_scales_key->append_text("C"); d_scales_key->append_text("C#"); d_scales_key->append_text("D");
            d_scales_key->append_text("D#");d_scales_key->append_text("E");  d_scales_key->append_text("F");
            d_scales_key->append_text("F#");d_scales_key->append_text("G");  d_scales_key->append_text("G#");
            d_scales_key->append_text("A"); d_scales_key->append_text("A#"); d_scales_key->append_text("B");
            d_scales_key->set_active(0);

            b->pack_start(*d_scales_key);
            d_scales_key->show();

            d_scales_mode->append_text("Chromatic");
            d_scales_mode->append_text("Major");
            d_scales_mode->append_text("Minor");
            d_scales_mode->append_text("Melodic minor 1");
            d_scales_mode->append_text("Melodic minor 2");
            d_scales_mode->append_text("Harmonic minor");
            d_scales_mode->append_text("Whole tone");
            d_scales_mode->append_text("Pentatonic 1");
            d_scales_mode->append_text("Pentatonic 2");
            d_scales_mode->set_active(0);

            b->pack_start(*d_scales_mode);
            d_scales_mode->show();

            d_scales_key->signal_changed().connect( sigc::mem_fun(*this, &VocProcGUI::scale_change) );
            d_scales_mode->signal_changed().connect( sigc::mem_fun(*this, &VocProcGUI::scale_change) );


            refBuilder->get_widget("sThreshold", s_threshold);
            s_threshold->set_range(0, 1);
            s_threshold->set_increments(0.01, 0.1);
            slot<void> slot_threshold = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_threshold),
                mem_fun(*s_threshold, &VScale::get_value));
            s_threshold->signal_value_changed().connect(slot_threshold);

            refBuilder->get_widget("sAttack", s_attack);
            s_attack->set_range(0, 1);
            s_attack->set_increments(0.01, 0.1);
            slot<void> slot_attack = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_attack),
                mem_fun(*s_attack, &VScale::get_value));
            s_attack->signal_value_changed().connect(slot_attack);

            refBuilder->get_widget("sTranspose", s_transpose);
            s_transpose->set_range(-12, 12);
            s_transpose->set_increments(1, 2);
            slot<void> slot_transpose = compose(bind<0>(mem_fun(*this, &VocProcGUI::write_control), p_transpose),
                mem_fun(*s_transpose, &VScale::get_value));
            s_transpose->signal_value_changed().connect(slot_transpose);

            // next line crashes ardour if enabled (ingen is OK with it)
            //scale_change();

            // -----------------

            add(*main_box);
        }
    }


    void scale_change(){

        int params[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        // intervals
        int chromatic[]={11, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        int major[]={7, 2, 2, 1, 2, 2, 2, 1};
        int minor[]={7, 2, 1, 2, 2, 1, 2, 2};
        int melodic_minor_1[]={7, 2, 1, 2, 2, 2, 2, 1};
        int melodic_minor_2[]={7, 2, 2, 1, 2, 2, 1, 2};
        int harmonic_minor[]={7, 2, 1, 2, 2, 1, 3, 1};
        int whole_tone[]={5, 2, 2, 2, 2, 2};
        int pentatonic_1[]={4, 2, 3, 2, 2};
        int pentatonic_2[]={4, 2, 2, 3, 2};


        int *mode;

        int idx=d_scales_key->get_active_row_number();

        switch(d_scales_mode->get_active_row_number()){
            case 0:
                mode=chromatic;
                break;
            case 1:
                mode=major;
                break;
            case 2:
                mode=minor;
                break;
            case 3:
                mode=melodic_minor_1;
                break;
            case 4:
                mode=melodic_minor_2;
                break;
            case 5:
                mode=harmonic_minor;
                break;
            case 6:
                mode=whole_tone;
                break;
            case 7:
                mode=pentatonic_1;
                break;
            case 8:
                mode=pentatonic_2;
                break;
            default:
                mode=chromatic;
        }

        params[idx]=1;
        for(int i=0;i<mode[0];i++){
            idx=(idx+mode[i+1])%12;
            params[idx]=1;
        }

        s_transpose->set_range(-mode[0]-1, mode[0]+1);
        if(s_transpose->get_value() < -mode[0]-1) s_transpose->set_value(-mode[0]-1);
        if(s_transpose->get_value() > mode[0]+1) s_transpose->set_value(mode[0]+1);

        // write controls
        for(int i=0;i<12;i++) write_control(p_c+i, params[i]);

    }

    void port_event(uint32_t port, uint32_t buffer_size, uint32_t format, const void* buffer) {
        float val=*static_cast<const float*>(buffer);
        switch(port){
            case p_pitch_factor:
                s_pitch_factor->set_value(val);
                break;

            case p_effect:
                s_effect->set_value(val);
                Label *l;
                refBuilder->get_widget("l_effect", l);
                if(val==0.0) l->set_text("effect off");
                else l->set_text("effect");
                break;

            case p_fc_voc_switch:
                if(val<0.5) c_for_voc->set_active(false);
                else c_for_voc->set_active(true);
                break;

#ifndef NO_VOCODER
            case p_fc_voc:
                d_for_voc->set_active((int)(val+0.5));
                break;
#endif

            case p_pitch_correction:
                if(val<0.5){
                    s_pitch_factor->set_sensitive(true);
                    c_autotune->set_active(false);
                }
                else{
                    s_pitch_factor->set_value(0);
                    s_pitch_factor->set_sensitive(false);
                    c_autotune->set_active(true);
                }
                break;

            case p_threshold:
                s_threshold->set_value(val);
                break;

            case p_attack:
                s_attack->set_value(val);
                break;

            case p_transpose:
                s_transpose->set_value(val);
                break;

            case p_offset:
                refBuilder->get_widget("pLeft", p_left);
                refBuilder->get_widget("pCenter", p_center);
                refBuilder->get_widget("pRight", p_right);
                if(val==-100){
                    p_left->set_fraction(1.0);
                    p_center->set_fraction(0.0);
                    p_right->set_fraction(0.0);
                }
                else if(val==100){
                    p_left->set_fraction(0.0);
                    p_center->set_fraction(0.0);
                    p_right->set_fraction(1.0);
                }
                else if(val<0){
                    p_left->set_fraction(-val/100.0);
                    p_center->set_fraction(1.0);
                    p_right->set_fraction(0.0);
                }
                else if(val>0){
                    p_left->set_fraction(0.0);
                    p_center->set_fraction(1.0);
                    p_right->set_fraction(val/100.0);
                }
                else if(val==0){
                    p_left->set_fraction(0.0);
                    p_center->set_fraction(1.0);
                    p_right->set_fraction(0.0);
                }
                

                break;
        }
        
    }


protected:
    
    Glib::RefPtr<Builder> refBuilder;

    VScale *s_pitch_factor;
    VScale *s_effect;

    VScale *s_attack;
    VScale *s_threshold;
    VScale *s_transpose;

    CheckButton *c_for_voc;
    CheckButton *c_autotune;

    ProgressBar *p_left;
    ProgressBar *p_right;
    ProgressBar *p_center;

    ComboBoxText *d_for_voc;
    ComboBoxText *d_scales_key;
    ComboBoxText *d_scales_mode;
};


static int _ = VocProcGUI::register_class("http://hyperglitch.com/dev/VocProc/gui");

