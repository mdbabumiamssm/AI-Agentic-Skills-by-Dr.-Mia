# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Sidebar with About section
import streamlit as st

def render_sidebar():
    st.sidebar.header("About")
    st.sidebar.info(
        """
        **AI Speech Trainer** helps users improve their public speaking skills through:\
        
        
        📽️ Video Analysis\
        
        🗣️ Voice Analysis\
                
        📊 Content Analysis & Feedback\
        

        - Upload your video to receive a detailed feedback.

        - Improve your public speaking skills with AI-powered analysis.
        
        - Get personalized suggestions to enhance your performance.
        """
    ) 
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
